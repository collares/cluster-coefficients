#ifndef CLUSTER_CARTPROD_BICLIQUE_INST_H
#define CLUSTER_CARTPROD_BICLIQUE_INST_H

#include "biclique.h"
#include "cartprod.h"

/*
  This instance computes cluster expansion coefficients for the cartesian
  product of t copies of K_{s,s}. For computing L_j (j >= 2), it does so by
  considering "embeddable" structures which are subsets of the cartesian
  product of 2*(j-1) copies of K_{j,j-1}.
*/
class cartprod_biclique_inst : public instance<cartprod_biclique_inst> {
private:
    int j_;
    GiNaC::symbol s_, t_;

public:
    using G = cartprod<biclique>;
    using V = cartprod<biclique>::V;
    G graph_;

    cartprod_biclique_inst(int j, GiNaC::symbol s, GiNaC::symbol t)
        : instance(j), j_(j), s_(s), t_(t), graph_(2*std::max(1, j-1)) { }

    V root() const {
        return V(std::vector<biclique_vertex>(graph_.num_coords(), 0));
    }

    struct data {
        int dist_to_embeddable;
        GiNaC::ex embeddings;
    };

    unsigned int
    active_union(const big_polymer<cartprod_biclique_inst>& bp) const {
        V rt = root();

        unsigned int ans = 0;
        for (int i = 0; i < bp.size(); i++)
            ans |= bp.vtx(i).active_coords(rt);
        return ans;
    }

    /*
      For us, a polymer is embeddable if:

      1) In each coordinate i, the set O_i (resp E_i) of odd (resp even)
      values seen in that coordinate is a prefix of the odd (resp even)
      positive integers.

      2) |O_i| + |E_i| >= |O_{i+1}| + |E_{i+1}| for every i. In particular,
      this implies that the set of active coordinates form a prefix.
    */

    data active_info(const big_polymer<cartprod_biclique_inst>& bp,
        bool compute_embeddings = false) const {
        GiNaC::ex emb = 1;

        int max_gap = 0;
        int max_size_when_filled = 0;
        int cur_run = 0;
        int dist = 0;
        int first_nonempty = -1;

        for (int j = graph_.num_coords() - 1; j >= 0; j--) {
            unsigned int coord_masks[2] = {0, 0};
            for (int i = 0; i < bp.size(); i++) {
                auto& v = bp.vtx(i);
                coord_masks[v.get(j).parity()] |= 1<<(v.get(j).data_ >> 1);
            }
            // the root is even but should not be counted
            assert(coord_masks[0] & 1);
            coord_masks[0] >>= 1;

            if (first_nonempty == -1) {
                if (!coord_masks[0] && !coord_masks[1])
                    continue;
                first_nonempty = j;
            }

            // Property 2 implies `size_when_filled` must appear in
            // non-decreasing order as we iterate backwards.
            int size_when_filled = std::bit_width(coord_masks[0])
                + std::bit_width(coord_masks[1]);
            int to_increasing = max_size_when_filled - size_when_filled;
            int cur_gap = size_when_filled
                - std::popcount(coord_masks[0])
                - std::popcount(coord_masks[1]);
            dist += cur_gap;

            if (to_increasing < 0) {
                max_size_when_filled = size_when_filled;
                max_gap = std::max(max_gap, cur_gap);
                cur_run = 1;
            } else {
                dist += to_increasing;
                max_gap = std::max(max_gap, cur_gap + to_increasing);
                cur_run++;
            }

            if (compute_embeddings) {
                emb *= (t_ - (first_nonempty - j)) / cur_run;
                emb *= GiNaC::binomial(s_ - 1, std::bit_width(coord_masks[0]));
                emb *= GiNaC::binomial(s_, std::bit_width(coord_masks[1]));
            }
        }

        assert(!compute_embeddings || dist == 0);
        return {std::max(max_gap, (dist+1)/2), emb};
    }

    std::string
    debug_data(const big_polymer<cartprod_biclique_inst>& bp) const {
        auto info = active_info(bp, true);
        std::stringstream ss;
        ss << "active mask = " << active_union(bp) << ", "
           << "number of unrooted embeddings = " << info.embeddings << ", "
           << "distance = " << info.dist_to_embeddable;
        return ss.str();
    }

    int
    dist_to_embeddable(const big_polymer<cartprod_biclique_inst>& bp) const {
        return active_info(bp).dist_to_embeddable;
    }

    bool is_embeddable(const big_polymer<cartprod_biclique_inst>& bp) const {
        return dist_to_embeddable(bp) == 0;
    }

    GiNaC::ex
    rooted_embeddings(const big_polymer<cartprod_biclique_inst>& bp) const {
        return pow(2*s_, t_) * active_info(bp, true).embeddings / 2;
    }

    inline bool is_valid(const subpolymer<cartprod_biclique_inst>& sp) const {
        return true;
    }

    /*
      This represents N_c(v). We could do it more efficiently by replacing
      coordinate c in vertex v by a sentinel (either 2j or 2j+1), but it would
      be a big abstraction leak and it doesn't seem worth it at the moment.

      TODO: This probably should be a part of cartprod.
    */
    struct neighbourhood {
    private:
        V vertex_;
        int changed_coord_;

    public:
        neighbourhood(V v, int changed_coord) :
            vertex_(v), changed_coord_(changed_coord) { }

        int new_parity() const {
            return !vertex_.get(changed_coord_).parity();
        }

        bool operator<(const neighbourhood& w) const {
            assert(vertex_.num_coords() == w.vertex_.num_coords());

            if (changed_coord_ != w.changed_coord_)
                return changed_coord_ < w.changed_coord_;

            if (new_parity() != w.new_parity())
                return new_parity() < w.new_parity();

            for (int i = 0; i < vertex_.num_coords(); i++)
                if (i != changed_coord_ && vertex_.get(i) != w.vertex_.get(i))
                    return vertex_.get(i) < w.vertex_.get(i);

            return false;
        }

        bool contains(const V& v) const {
            if (new_parity() != v.get(changed_coord_).parity())
                return false;

            for (int i = 0; i < vertex_.num_coords(); i++)
                if (i != changed_coord_ && vertex_.get(i) != v.get(i))
                    return false;

            return true;
        }

        /*
          The first coordinate of the return value is true if the intersection
          is a single vertex, in which case the second coordinate is the
          intersection. The second coordinate is to be ignored if the first one
          equals false. Would be nice to have algebraic data types so we could
          return the intersetion properly, but this will do.
        */
        std::pair<bool, V> intersect(const neighbourhood& n) const {
            if (changed_coord_ == n.changed_coord_ ||
                n.vertex_.get(changed_coord_).parity() != new_parity())
                return {false, vertex_};

            V candidate = vertex_;
            candidate = candidate.set(changed_coord_,
                                      n.vertex_.get(changed_coord_));
            assert(contains(candidate));

            return {n.contains(candidate), candidate};
        }
    };

    /*
      Recall that if the polymer is P, the weight of the corresponding polymer
      is lambda^|P| / (1+lambda)^{|N(P)|}. Therefore, the weight function is
      essentially a function that computes |N(P)|. In our case, we shall:

      1) Compute a set S of all w such that $N_{c_1}(v_1) \cap N_{c_2}(v_2) =
      {w}$ for two pairs (v1, c1) and (v2, c2). c1 and c2 must be active. This
      set has size O(j^3), since each element must differ from one of the j
      elements of P in one coordinate. The size of S will be added to the final
      answer.

      2) Compute all the family of sets T of the form N_{c}(v) with c active.
      We do so to avoid double counting, since if v_1 and v_2 differ only in
      coordinate c and have the same parity there, then $N_{c}(v_1) =
      N_c(v_2)$, but we do not want to count those twice. Our goal is to say
      each of these sets contribute s to the neighbourhood, but for each of
      these sets we must subtract its intersection with the set S from step 1.

      We don't try to be clever in step 2: for each set T, we naively intersect
      it with S.

      3) Add the contribution from the inactive coordinates. If there are "a"
      active coordinates, that is just num_vtxs * s_ * (t_ - a).
    */
    inline GiNaC::ex weight(const subpolymer<cartprod_biclique_inst>& sp,
                            GiNaC::symbol lambda) const {
        int num_vtxs = std::popcount(sp.mask_);
        int a = std::popcount(active_union(sp.bp_));
        GiNaC::ex ans;

        // step 1
        std::set<V> S;
        for (int it = sp.mask_; it; it &= it-1) {
            V v1 = sp.bp_.vtx(__builtin_ctz(it));
            for (int c1 = 0; c1 < a; c1++) {
                neighbourhood n1(v1, c1);
                for (int it2 = it; it2; it2 &= it2-1) {
                    V v2 = sp.bp_.vtx(__builtin_ctz(it2));
                    for (int c2 = 0; c2 < a; c2++) {
                        auto inter = n1.intersect(neighbourhood(v2, c2));
                        if (inter.first)
                            S.insert(inter.second);
                    }
                }
            }
        }
        ans += S.size();

        // step 2
        std::set<neighbourhood> T;
        for (int it = sp.mask_; it; it &= it-1) {
            V v = sp.bp_.vtx(__builtin_ctz(it));

            for (int c = 0; c < a; c++) {
                neighbourhood n(v, c);
                if (T.find(n) != T.end())
                    continue;

                int overcount = 0;
                for (auto w : S)
                    if (n.contains(w))
                        overcount++;

                ans += s_ - overcount;
                T.insert(n);
            }
        }

        // step 3
        ans += num_vtxs * s_ * (t_ - a);
        return pow(lambda, num_vtxs) / pow(1 + lambda, ans);
    }
};

#endif
