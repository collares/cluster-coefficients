/*
  This file contains the generic code necessary to compute cluster expansion
  coefficients. Specific graph implementations are in their own separate files.
*/

#ifndef CLUSTER_CLUSTER_H
#define CLUSTER_CLUSTER_H

#include <cassert>
#include <ginac/ginac.h>
#include <iostream>
#include <queue>
#include <set>
#include <vector>

template<typename I> class big_polymer;
template<typename I> std::ostream &operator<<(std::ostream &,
                                              const big_polymer<I> &);

/*
  A big polymer is simply a 2-linked ordered set of vertices containing the
  root.

  FIXME: Due to design decisions which are now obsolete, big_polymer does not
  support adding/deleting vertices, which causes a lot of copying when
  generating big polymers in the first place, since dist_to_embeddable and
  is_embeddable take a big_polymer argument and not a "partial polymer" one.
*/

template <typename I> class big_polymer {
private:
    std::vector<typename I::V> data;

public:
    big_polymer(const std::vector<typename I::V>& shv) : data(shv) { }

    const typename I::V vtx(int idx) const {
        return data[idx];
    }

    size_t size() const {
        return data.size();
    }

    friend std::ostream& operator<< <I>(std::ostream& o,
                                        const big_polymer<I>& bp);
};

template <typename I>
std::ostream& operator<<(std::ostream& o, const big_polymer<I>& bp) {
    bool first = true;
    o << "{";
    for (auto vtx : bp.data) {
        if (!first)
            o << ", ";
        o << vtx;
        first = false;
    }
    return o << "}";
}

/*
  Encodes a subset of a big polymer. Currently implemented as a mask.
*/
template<typename I>
class subpolymer {
public:
    const big_polymer<I>& bp_;
    unsigned int mask_;

    subpolymer(const big_polymer<I>& bp, unsigned int mask = 0)
        : bp_(bp), mask_(mask) {
    }

    // Two polymers are incompatible if their union is 2-linked.
    bool is_incompatible_with(const subpolymer<I>& other) const {
        for (int it1 = mask_; it1; it1 &= it1-1) {
            typename I::V v = bp_.vtx(__builtin_ctz(it1));
            for (int it2 = other.mask_; it2; it2 &= it2-1) {
                typename I::V w = other.bp_.vtx(__builtin_ctz(it2));
                if (v.dist(w) <= 2)
                    return true;
            }
        }

        return false;
    }

    size_t size() const {
        return std::popcount(mask_);
    }

    std::vector<typename I::V> vertices() const {
        std::vector<typename I::V> ans;
        for (int it = mask_; it; it &= it-1)
            ans.push_back(bp_.vtx(__builtin_ctz(it)));
        return ans;
    }

    bool is_two_linked() const {
        assert(mask_);

        std::queue<int> q;
        int visited_mask;
        int start = __builtin_ctz(mask_);

        q.push(start);
        visited_mask = 1<<start;
        int comp_sz = 1;

        while (!q.empty()) {
            int idx = q.front();
            q.pop();

            for (int it = mask_ & ~visited_mask; it; it &= it-1) {
                int nidx = __builtin_ctz(it);
                if (bp_.vtx(idx).dist(bp_.vtx(nidx)) == 2) {
                    q.push(nidx);
                    visited_mask |= 1<<nidx;
                    comp_sz++;
                }
            }
        }

        return comp_sz == size();
    }
};

template<typename I>
std::ostream& operator<<(std::ostream& o, const subpolymer<I>& sp) {
    bool first = true;
    o << "{";
    for (int it = sp.mask_; it; it &= it-1) {
        int idx = __builtin_ctz(it);
        if (!first)
            o << ", ";
        o << sp.bp_.vtx(idx);
        first = false;
    }
    return o << "}";
}

template<typename I> class cluster;
template<typename I> std::ostream &operator<<(std::ostream &,
                                              const cluster<I> &);

/*
  A cluster is an ordered collection of subpolymers of a big polymer. Each
   vertex of the big polymer must be in at least one subpolymer.
*/
template<typename I>
class cluster {
private:
    const big_polymer<I>& bp_;
    size_t total_sz_;
    std::stack<unsigned int> union_mask_;

public:
    std::vector<subpolymer<I>> elems_;

    cluster(const big_polymer<I>& bp) : bp_(bp), total_sz_(0) {
        union_mask_.push(0);
    }

    void push_subpolymer(const subpolymer<I>& sp) {
        assert(&bp_ == &sp.bp_);
        total_sz_ += sp.size();
        elems_.push_back(sp);
        union_mask_.push(union_mask_.top() | sp.mask_);
    }

    void pop_subpolymer() {
        total_sz_ -= elems_.back().size();
        elems_.pop_back();
        union_mask_.pop();
    }

    bool covers_big_polymer() {
        return union_mask_.top()+1 == 1<<bp_.size();
    }

    unsigned int uncovered() {
        return ((1<<bp_.size())-1) & ~(union_mask_.top());
    }

    size_t total_size() {
        return total_sz_;
    }

    size_t num_polymers() const {
        return elems_.size();
    }

    /*
      For a graph H represented by an adjacency matrix (in bitmask format), the
      function below computes (-1)^{v(H)-1} T_H(1, 0) / v(H)! by using the
      deletion-contraction formula. No loops or multiple edges are ever
      present.

      TODO: Produce a canonical label of the input graph. This is ideal, since
      there are only 208 graphs with at most 6 vertices up to isomorphism, and
      the number of polymers (vertices of the incompatibility graph) is not
      big in practice.

      TODO: This should not be here, but in a global cache.
    */
    GiNaC::ex
    delete_contract(int n,
                    std::vector<unsigned int>& adj_masks) const {
        // `cache_msk` needs to have enough bits to represent all (n choose 2)
        // edges, so we assert n <= 8. We can bump it to n <= 11 if `cache_msk`
        // is 64-bit.
        assert(n >= 1 && n <= 8);
        if (n == 1)
            return 1;

        unsigned int cache_msk = 0;
        for (int i = 0; i < n; i++) {
            cache_msk <<= n-1-i;
            cache_msk |= adj_masks[i]>>(i+1);
        }

        if (std::popcount(cache_msk) + 1 < n)
            return 0;

        static std::map<std::pair<int, unsigned int>, GiNaC::ex> cache;
        auto it = cache.find(std::make_pair(n, cache_msk));
        if (it != cache.end())
            return it->second;

        int ev = -1, ew = -1;
        {
            int pretime = 0;
            int num[n], low[n];
            memset(num, -1, sizeof num);

            std::function<void(int, int)> dfs = [&](int v, int parent) {
                num[v] = low[v] = pretime++;
                for (unsigned int it = adj_masks[v]; it; it &= it-1) {
                    int w = __builtin_ctz(it);
                    if (w == parent)
                        continue;

                    if (num[w] == -1) {
                        dfs(w, v);
                        low[v] = std::min(low[v], low[w]);
                    } else
                        low[v] = std::min(low[v], num[w]);

                    if (low[w] > num[v]) {
                        ev = v;
                        ew = w;
                    }
                }
            };

            dfs(0, -1);
            // The graph is disconnected. This never happens in recursive calls.
            if (pretime != n)
                return cache[std::make_pair(n, cache_msk)] = 0;
        }

        if (ew < ev)
            std::swap(ev, ew);

        GiNaC::ex ans;
        if (ev == -1) {
            // No edge is a bridge, so we need to delete as well as contract.
            for (int v = 0; v < n; v++)
                if (adj_masks[v]) {
                    ev = v;
                    ew = __builtin_ctz(adj_masks[v]);
                    break;
                }

            adj_masks[ev] &= ~(1U<<ew);
            adj_masks[ew] &= ~(1U<<ev);
            ans = delete_contract(n, adj_masks);
            adj_masks[ev] |= 1U<<ew;
            adj_masks[ew] |= 1U<<ev;
        }

        // We will contract the edge {ev, ew}. `ew` will be deleted. We use
        // that ev is not n-1.
        {
            std::vector<unsigned int> cont_masks = adj_masks;
            cont_masks[ev] &= ~(1U<<ew);
            cont_masks[ew] &= ~(1U<<ev);

            for (unsigned int it = cont_masks[ew]; it; it &= it-1) {
                int w2 = __builtin_ctz(it);
                cont_masks[ev] |= 1U<<w2;
                cont_masks[w2] |= 1U<<ev;
                cont_masks[w2] &= ~(1U<<ew);
            }
            cont_masks[ew] = 0;

            if (ew != n-1) {
                for (unsigned int it = cont_masks[n-1]; it; it &= it-1) {
                    int w2 = __builtin_ctz(it);
                    cont_masks[w2] |= 1U<<ew;
                    cont_masks[w2] &= ~(1U<<(n-1));
                }
                std::swap(cont_masks[ew], cont_masks[n-1]);
            }

            ans -= delete_contract(n-1, cont_masks)/n;
        }

        return cache[std::make_pair(n, cache_msk)] = ans;
    }

    /*
       The Ursell function of a graph is the sum over all subsets of edges A
       inducing a connected graph of (-1)^{|A|}, normalized by 1/v(H)!.
       Previous tests showed that even iterating over all subsets of edges,
       using the straightforward definition, was fast enough as long as we
       cached answers. But, since it equals (-1)^{v(H)-1} T_H(1, 0) / v(H)!,
       where T_H is the Tutte polynomial of the graph, we might as well compute
       the desired quantity using deletion contraction.
    */
    GiNaC::ex compute_ursell() const {
        int k = num_polymers();
        if (k == 1)
            return 1;

        std::vector<unsigned int> adj_masks(k);
        for (int i = 0; i < k; i++)
            for (int j = i+1; j < k; j++)
                if (elems_[i].is_incompatible_with(elems_[j])) {
                    adj_masks[i] |= 1U<<j;
                    adj_masks[j] |= 1U<<i;
                }

        return delete_contract(k, adj_masks);
    }

    /*
      A naive function which just iterates over all subsets of edges. For
      verification purposes only, but for small graphs it is quite fast.
      It also asserts that deletion-contraction computes the same answer.
    */
    GiNaC::ex compute_ursell_brute_force() const {
        assert(num_polymers());
        if (num_polymers() == 1)
            return 1;

        static std::map<std::pair<int, unsigned int>, GiNaC::ex> cache;

        unsigned int cache_msk = 0;
        // Should use larger type for mask for more polymers.
        assert(num_polymers() <= 8);
        for (int i = 0, msk_pos = 0; i < num_polymers(); i++)
            for (int j = i+1; j < num_polymers(); j++, msk_pos++)
                if (elems_[i].is_incompatible_with(elems_[j]))
                    cache_msk |= 1U<<msk_pos;

        if (std::popcount(cache_msk) + 1 < num_polymers()) {
            assert(compute_ursell() == 0);
            return 0;
        }

        auto it = cache.find(std::make_pair(num_polymers(), cache_msk));
        if (it != cache.end())
            return it->second;

        GiNaC::ex ans;
        std::vector<std::pair<int, int>> edges;

        for (int i = 0, msk_pos = 0; i < num_polymers(); i++)
            for (int j = i+1; j < num_polymers(); j++, msk_pos++)
                if (cache_msk & (1<<msk_pos))
                    edges.emplace_back(i, j);

        for (unsigned int mask = (1<<edges.size())-1; mask; mask--) {
            if (std::popcount(mask) + 1 < num_polymers())
                continue;

            std::vector<std::vector<int>> adj(num_polymers());

            for (int it = mask; it; it &= it-1) {
                int i = __builtin_ctz(it);
                adj[edges[i].first].push_back(edges[i].second);
                adj[edges[i].second].push_back(edges[i].first);
            }

            std::queue<int> q;
            std::vector<bool> visited(num_polymers());
            q.push(0);
            visited[0] = true;

            int comp_sz = 1;
            while (!q.empty()) {
                int idx = q.front();
                q.pop();

                for (auto nidx : adj[idx])
                    if (!visited[nidx]) {
                        visited[nidx] = true;
                        q.push(nidx);
                        comp_sz++;
                    }
            }

            if (comp_sz == num_polymers())
                ans += std::popcount(mask) % 2 ? -1 : 1;
            else if (mask == (1<<edges.size())-1)
                break;
        }

        assert(compute_ursell() == ans / GiNaC::factorial(num_polymers()));
        return cache[std::make_pair(num_polymers(), cache_msk)] =
            ans / GiNaC::factorial(num_polymers());
    }

    friend std::ostream& operator<< <I>(std::ostream& o, const cluster<I>& c);
};

template<typename I>
std::ostream& operator<<(std::ostream& o, const cluster<I>& c) {
    bool first = true;
    o << "(";
    for (auto sp : c.elems_) {
        if (!first)
            o << ", ";
        o << sp;
        first = false;
    }
    return o << ")";
}

/*
  A design goal of this program is that computing different coefficients for
  the same graph uses the same generated code to maximize flexibility, but
  different graphs use different generated code (via template monomorphisation)
  for efficiency. For this reason, we're using the CRTP for static dispatch and
  inlining.

  TODO: This entails a bit of boilerplate. Check if it really speeds things up,
  or even if it is needed. When C++23 is more widespread, we can use "deducing
  this".
*/
template<typename I> class instance {
private:
    GiNaC::ex ans_;

    unsigned long long processed_tuples_;
    unsigned long long processed_multisets_;
    unsigned long long processed_big_polymers_;

    int j_; // this instance computes L_j

public:
    instance(int j) : j_(j) { }

    /*
      In the hypercube case, this is true if and only if the active coordinates
      form a prefix. This requires knowing the root.
    */
    inline bool is_embeddable(const big_polymer<I>& bp) const {
        return static_cast<const I*>(this)->is_embeddable(bp);
    }

    /*
      A lower bound on the number of vertices (each taken from the second
      neighbourhood of a previous vertex of bp) needed to add to bp to make it
      embeddable, used as a heuristic in the big polymer generation process.
      For now, non-embeddable polymers can return 0 (we use the is_embeddable
      function to be sure they are embeddable).
    */
    inline int dist_to_embeddable(const big_polymer<I>& bp) const {
        return static_cast<const I*>(this)->dist_to_embeddable(bp);
    }

    /*
      In the hypercube case, the factor is \binom{d}{a} / j whenever the big
      polymer has j vertices and active coordinate set [a].

      This should be zero when !is_embeddable(bp), and ideally it should be
      nonzero otherwise.
    */
    inline GiNaC::ex rooted_embeddings(const big_polymer<I>& bp) const {
        return static_cast<const I*>(this)->rooted_embeddings(bp);
    }

    inline std::string debug_data(const big_polymer<I>& bp) const {
        return static_cast<const I*>(this)->debug_data(bp);
    }

    /*
      It might be the case that not every subset of vertices of a big polymer
      is a valid polymer. In the antichains case, for example, polymers have to
      live entirely in one layer, even though the big polymer spans two layers.
    */
    inline bool is_valid(const subpolymer<I>& sp) const {
        return static_cast<const I*>(this)->is_valid(sp);
    }

    /*
      We've computed the polymers with regards to an induced subgraph G of the
      d-dimensional structure (d is symbolic) which (up to symmetry) is large
      enough to contain every vertex we need to consider to compute the desired
      coefficients. This function computes the weight of the polymer as a
      function of d.
    */
    inline GiNaC::ex weight(const subpolymer<I>& sp,
                            GiNaC::symbol lambda) const {
        return static_cast<const I*>(this)->weight(sp, lambda);
    }

    void process_big_polymer(const big_polymer<I>& bp,
                             GiNaC::symbol lambda) {
        // std::cerr << "* Processing big polymer " << bp
        //           << ", root = " << static_cast<const I*>(this)->root()
        //           << ", " << debug_data(bp) << std::endl;
        std::vector<subpolymer<I>> sps;
        std::vector<GiNaC::ex> weights;

        /*
          Generate all possible 2-linked subpolymers (polymers are always
          non-empty).

          TODO: The clever thing to do here would be to generate the
          incompatibility graph for all subpolymers here (instead of
          computing for the chosen subpolymers when computing the Ursell
          function), and use this to only generate clusters for which the
          incompatibility graph is connected.
        */
        for (int mask = 1; mask < (1 << bp.size()); mask++)
            if (mask & (mask-1)) {
                subpolymer<I> sp(bp, mask);
                if (sp.is_two_linked() && is_valid(sp)) {
                    sps.push_back(sp);
                    weights.push_back(weight(sp, lambda));
                }
            }

        // Leaving the singletons for last ensures the backtacking never
        // reaches a dead end.
        for (int mask = 1; mask < (1 << bp.size()); mask <<= 1) {
            subpolymer<I> sp(bp, mask);
            if (is_valid(sp)) {
                sps.push_back(sp);
                weights.push_back(weight(sp, lambda));
            }
        }

        /*
          Now we want to find all ways of picking clusters, which are
          tuples of elements of sps whose union is bp (possibly with
          repeats) and whose sum of sizes is j_. Since the weight of a
          cluster does not depend on the ordering of its polymers, we may
          do that by finding multisets of polymers with the above property
          then multiplying by the appropriate number of permutations.
        */

        GiNaC::ex bp_contrib = 0;
        cluster<I> cl(bp); // start with the empty cluster, then backtrack
        std::vector<int> ids;

        std::function<void(int, int, int)> backtrack =
            [&](int last_idx, int cur_run, int num_tuples) {
                if (cl.total_size() == j_ && cl.covers_big_polymer()) {
                    // We could have a cl.weight() function, but caching the
                    // weights of subpolymers is significantly more efficient
                    // in some cases and this is the simplest way to do so.
                    GiNaC::ex cur = cl.compute_ursell();
                    for (auto id : ids)
                        cur *= weights[id];

                    bp_contrib += num_tuples * cur;
                    processed_tuples_ += num_tuples;
                    processed_multisets_++;
                    // std::cerr << "** Cluster " << cl << " has weight "
                    //            << cur << " and permutations of it occur "
                    //            << num_tuples << " times " << std::endl;
                    return;
                }

                if (cl.total_size() + std::popcount(cl.uncovered()) > j_)
                    return;

                // It would be good to improve this part to ensure polynomial
                // delay.
                for (int i = last_idx; i < sps.size(); i++) {
                    int next_run = last_idx == i ? cur_run + 1 : 1;
                    cl.push_subpolymer(sps[i]);
                    ids.push_back(i);
                    backtrack(i,
                              next_run,
                              cl.num_polymers() * num_tuples / next_run);
                    ids.pop_back();
                    cl.pop_subpolymer();

                    // Subpolymers of size 1 are tried last. If we skip the
                    // last chance of covering a vertex, no need to continue.
                    if (sps[i].size() == 1 && (cl.uncovered() & sps[i].mask_))
                        break;
                }
            };

        backtrack(0, 0, 1);
        ans_ += bp_contrib * rooted_embeddings(bp) / bp.size();
    }

    /*
      Generates the set of embeddable big polymers and processes them.

      A "big polymer" is a 2-linked set which lives on the defect side and
      which will be the vertex set of a cluster. By definition, we also require
      all big polymers to be _embeddable_, that is, that its set of "active
      coordinates" (i.e. coordinates which differ from the fixed root) satisfy
      a simple structure which allows us to map them to the supergraph.
      Concretely, this means that the set of i->(1-i) active coordinates is a
      prefix of the originally-i coordinates.

      j_ is the size of the cluster. Most of this function generates big
      polymers of tize at most j_. It then delegates the task of processing
      them to the function `process_big_polymer` above.
    */
    GiNaC::ex compute(GiNaC::symbol lambda) {
        assert(j_ < 31);

        ans_ = 0;
        processed_big_polymers_ = 0;
        processed_tuples_ = 0;
        processed_multisets_ = 0;

        // Invariant: seen_set and seen contain the same elements.
        std::vector<typename I::V> seen;
        std::set<typename I::V> seen_set;

        typename I::G graph = static_cast<const I*>(this)->graph_;
        typename I::V rt = static_cast<const I*>(this)->root();
        seen_set.insert(rt);
        seen.push_back(rt);

        // cur is always 2-linked, but may not be embeddable. As the
        // backtracking progresses, it runs over all 2-linked subsets of size k
        // containing rt whose dist_to_embeddable is at most <= j_ - k (for 1
        // <= k <= j_).
        std::vector<typename I::V> cur;

        std::function<void(int)> backtrack = [&](int fringe_start) {
            // Invariant: Every element of seen[0], ..., seen[fringe_start - 1]
            // has been explored (that is, it is either in cur or will never be
            // considered in subcalls). The remaining elements in `seen` (the
            // "fringe") are still to be processed. The fringe acts as a queue,
            // so we are trying possibilities breadth-first, despite what the
            // recursion might suggest.
            //
            // Every element of `seen` is a neighbour (in the second power) of
            // some element of cur. In particular, seen always has size
            // polynomial in j_. The recursion stack is cur.size() layers deep
            // at any given point.
            int seen_sz = seen.size();
            for (int i = fringe_start; i < seen_sz; i++) {
                auto& v = seen[i];
                cur.push_back(v);

                big_polymer<I> bp(cur);
                if (cur.size() <= j_ && is_embeddable(bp)) {
                    process_big_polymer(bp, lambda);
                    processed_big_polymers_++;
                    if ((processed_big_polymers_ & (1ULL<<18)-1) == 0)
                        std::cerr << "Processed " << processed_big_polymers_
                                  << " big polymers, last one was "
                                  << bp << std::endl;
                }

                if (cur.size() < j_ &&
                    dist_to_embeddable(bp) <= j_ - cur.size()) {
                    for (auto& w : graph.second_neighbourhood(v))
                        if (seen_set.find(w) == seen_set.end()) {
                            seen.push_back(w);
                            seen_set.insert(w);
                        }

                    backtrack(i+1);

                    // Remove the elements we added above.
                    while (seen.size() > seen_sz) {
                        seen_set.erase(seen.back());
                        seen.pop_back();
                    }
                }

                cur.pop_back();
            }
        };

        backtrack(0);
        std::cerr << "Done, processed "
                  << processed_big_polymers_ << " big polymers and "
                  << processed_tuples_ << " clusters"
                  << " (" << processed_multisets_ << " multisets)"
                  << std::endl << std::endl;
        return ans_;
    }
};

#endif
