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
*/

template <typename I> class big_polymer {
private:
    // Annoyingly, during construction it is really convenient to have this as
    // a set, but we want to use the ordering later.
    std::vector<typename I::V> data;

public:
    big_polymer(const std::set<typename I::V>& shv) {
        data = std::vector(shv.begin(), shv.end());
    }

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
    big_polymer<I>& bp_;
    unsigned int mask_;

    subpolymer(big_polymer<I>& bp, unsigned int mask = 0)
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

    size_t total_size() {
        return total_sz_;
    }

    size_t num_polymers() const {
        return elems_.size();
    }

    /* The Ursell function of a graph is the sum over all subsets of edges A
       inducing a connected graph of (-1)^{|A|}, normalized by 1/v(H)!. This
       computes the Ursell function of the incompatibility graph of the
       cluster.

       The above implementation is very naive, but with caching it does not
       seem to be the bottleneck. Two approaches to improve it would be:

       1) Produce a canonical label of the input graph. This is ideal, since
       there are only 208 graphs with at most 6 vertices up to isomorphism, and
       the number of polymers (vertices of the incompatibility graph) is not
       big in practice.

       2) Use a smarter algorithm. Deletion/contraction (perhaps splitting into
       biconnected components) is already a theoretical improvement, but
       Jenssen--Perkins note that Ursell function of a graph is, up to a
       multiplicative factor of ((-1)^(1-v(H)))/v(H)!), an evaluation of the
       Tutte polynomial at (1, 0). Tutte polynomials may be evaluated in
       vertex-exponential time by an algorithm of Björklund, Husfeldt, Kaski
       and Koivisto. Both approaches are unlikely to improve the running time
       in practice, since we are interested in very small graphs (up to 6
       vertices).
    */
    GiNaC::ex compute_ursell() const {
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

        if (std::popcount(cache_msk) + 1 < num_polymers())
            return 0;

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
public:
    int j_; // this instance computes L_j

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

    // TODO: This is here because it relies on the specific weight function.
    // Should it really be here instead of on cluster, though?
    inline GiNaC::ex compute_weight(const cluster<I>& cl,
                                    GiNaC::symbol lambda) const {
        GiNaC::ex ans = cl.compute_ursell();
        for (auto sp : cl.elems_)
            ans *= weight(sp, lambda);
        return ans;
    }

    /*
      Generates the set of embeddable big polymers.

      A "big polymer" is a 2-linked set which lives on the defect side and
      which will be the vertex set of a cluster. By definition, we also require
      all big polymers to be _embeddable_, that is, that its set of "active
      coordinates" (i.e. coordinates which differ from the fixed root) satisfy
      a simple structure which allows us to map them to the supergraph.
      Concretely, this means that the set of i->(1-i) active coordinates is a
      prefix of the originally-i coordinates.

      TODO: Unfortunately the procedure described in Jenssen--Perkins can
      generate any big polymer multiple times, so we have to use a lot of
      memory for deduplication. There must be a fast way of generating each
      embeddable set once without using too much memory.
    */
    std::vector<big_polymer<I>>
    generate_embeddable_big_polymers(int max_elements) const {
        typename I::G graph = static_cast<const I*>(this)->graph_;
        using partial_polymer = std::set<typename I::V>;

        std::vector<big_polymer<I>> ret;

        partial_polymer start;
        start.insert(static_cast<const I*>(this)->root());
        if (is_embeddable(start))
            ret.push_back(start);

        std::set<partial_polymer> Lprev, Lcur;
        Lprev.insert(start);

        for (int k = 1; k < max_elements; k++) {
            Lcur.clear();

            for (const partial_polymer &prev_set : Lprev)
                for (auto v : prev_set) {
                    auto new_vtxs = graph.second_neighbourhood(v);
                    for (auto w : new_vtxs) {
                        partial_polymer new_set = prev_set;

                        if (new_set.find(w) == new_set.end()) {
                            new_set.insert(w);
                            assert(new_set.size() == k + 1);

                            big_polymer<I> bp(new_set);
                            if (Lcur.find(new_set) == Lcur.end() &&
                                dist_to_embeddable(bp) <= max_elements - (k+1))
                                Lcur.insert(new_set);
                        }
                    }
                }

            for (const partial_polymer &pp : Lcur) {
                big_polymer<I> bp(pp);
                if (is_embeddable(bp))
                    ret.push_back(bp);
            }

            swap(Lprev, Lcur);
        }

        return ret;
    }

    GiNaC::ex compute(GiNaC::symbol lambda) {
        GiNaC::ex ans;

        assert(j_ < 31);
        std::cerr << "Generating big polymers..." << std::endl;
        std::vector<big_polymer<I>> L =
            generate_embeddable_big_polymers(j_);

        unsigned long long processed_tuples = 0;
        unsigned long long processed_multisets = 0;

        std::cerr << "About to process " << L.size() << " big polymers"
                  << std::endl;
        for (auto &bp : L) {
            // std::cerr << "* Processing big polymer " << bp
            //           << ", root = " << static_cast<const I*>(this)->root()
            //           << ", " << debug_data(bp) << std::endl;

            std::vector<subpolymer<I>> sps;

            /*
              Generate all possible 2-linked subpolymers (polymers are always
              non-empty).

              TODO: The clever thing to do here would be to generate the
              incompatibility graph for all subpolymers here (instead of
              computing for the chosen subpolymers when computing the Ursell
              function), and use this to only generate clusters for which the
              incompatibility graph is connected.
            */
            for (int mask = 1; mask < (1 << bp.size()); mask++) {
                subpolymer sp(bp, mask);
                if (sp.is_two_linked() && is_valid(sp))
                    sps.push_back(sp);
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

            std::function<void(int, int, int)> backtrack =
                [&](int last_idx, int cur_run, int num_tuples) {
                    if (cl.total_size() == j_ && cl.covers_big_polymer()) {
                        GiNaC::ex cur = compute_weight(cl, lambda);
                        bp_contrib += num_tuples * cur;
                        processed_tuples += num_tuples;
                        processed_multisets++;
                        // std::cerr << "** Cluster " << cl << " has weight "
                        //            << cur << " and permutations of it occur "
                        //            << num_tuples << " times " << std::endl;
                    }

                    if (cl.total_size() >= j_)
                        return;

                    for (int i = last_idx; i < sps.size(); i++) {
                        int next_run = last_idx == i ? cur_run + 1 : 1;
                        cl.push_subpolymer(sps[i]);
                        backtrack(i,
                                  next_run,
                                  cl.num_polymers() * num_tuples / next_run);
                        cl.pop_subpolymer();
                    }
                };

            backtrack(0, 0, 1);
            ans += bp_contrib * rooted_embeddings(bp) / bp.size();
        }

        std::cerr << "Done, processed " << processed_tuples << " clusters"
                  << " (" << processed_multisets << " multisets)"
                  << std::endl << std::endl;
        return ans;
    }
};

#endif
