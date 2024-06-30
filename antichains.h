/*
  The graph has size n = 2d or n = 2d + 1, depending on whether odd = true or
  not.

  In the even case, the middle layer M is unique. In the odd case we're we will
  let M be the d-th layer. The (d+1)-th layer will be called U, and the
  (d-1)-th layer will be called L. In the odd case, L will be smaller than U.
  Together, L and U are the defect side. Note that since there are two choices
  for M in the odd case, the answer has to be multiplied by 2 at the end, but
  this is out of the scope of this program since the correction is not in the
  exponent.

  Polymers are still 2-linked subsets, and clusters are still tuples of
  polymers whose "union" is 2-linked. However, a "valid" polymer must be
  entirely contained in one of the layers. The same restriction does not apply
  to clusters), so to make sure every cluster whose underlying set has k vertex
  is counted k times, we need to fix a vertex in L, and later fix a vertex in
  U.
*/

#ifndef CLUSTER_ANTICHAINS_H
#define CLUSTER_ANTICHAINS_H

#include <bit>
#include <iostream>
#include <ginac/ginac.h>
#include <sstream>
#include "cluster.h"
#include "hypercube_vertex.h"

/*
  For d = floor(n/2), this represents layers d-1, d and d+1 of a n-dimensional
  hypercube. sz is always even in the applications, but this is not required.
*/
class middle_layers {
private:
    int sz_;

public:
    middle_layers(int sz) : sz_(sz) {
    }

    std::vector<hypercube_vertex> second_neighbourhood(hypercube_vertex v) {
        std::vector<hypercube_vertex> ans;
        for (int i = 0; i < sz_; i++)
            for (int j = i+1; j < sz_; j++) {
                auto w = v.flip(i).flip(j);
                if (abs(w.layer() - sz_/2) <= 1)
                    ans.push_back(w);
            }
        return ans;
    }
};

class antichains_inst : public instance<antichains_inst> {
private:
    int middle_ones_;
    int root_ones_;

    // symbolic variable corresponding to ambient hypercube dimension
    GiNaC::ex n_;

    // a vertex v will be mapped to v.layer() + mapping_delta in the
    // n-dimensional object.
    GiNaC::ex mapping_delta_;

    struct data {
        int dist_to_embeddable;
        int changed_zeros, changed_ones;
        unsigned int special_coords;
    };

    inline data prefix_info(const big_polymer<antichains_inst>& bp) const {
        unsigned int aunion = 0;
        int dist = 0;
        for (int i = 0; i < bp.size(); i++)
            aunion |= bp.vtx(i).active_coords(root_);

        unsigned int ones_to_zeros = ((1<<root_ones_) - 1) & aunion;
        dist += std::popcount(~ones_to_zeros) - std::countl_zero(ones_to_zeros);

        unsigned int zeros_to_ones = aunion >> root_ones_;
        dist += std::popcount(~zeros_to_ones) - std::countl_zero(zeros_to_ones);
        zeros_to_ones <<= root_ones_;

        return data{(dist+1)/2,
                    std::popcount(zeros_to_ones),
                    std::popcount(ones_to_zeros),
                    ones_to_zeros | zeros_to_ones};
    }

public:
    using G = middle_layers;
    using V = hypercube_vertex;
    middle_layers graph_;
    V root_;

    /*
      We can take our subgraph to always have an even dimension. j-1 "moves"
      starting from the lower layer require at least 2*(j-1) zeros (consider a
      star) and 2*(j-2) ones (move to the upper layer and consider a star).
      Therefore, the minimum dimension is 4*(j-1). We will use 4*j to not have
      to special case j==1. required. The middle layers (2*j-1, 2*j, 2*j+1)
      will be mapped to (L, M, U).
    */
    antichains_inst(int delta_middle,
                    int j,
                    GiNaC::ex n,
                    GiNaC::ex middle_layer)
        : instance(j),
          middle_ones_(2*j),
          root_ones_(middle_ones_ + delta_middle),
          n_(n),
          mapping_delta_(middle_layer - middle_ones_),
          graph_(middle_layers(2*middle_ones_)),
          root_(V((1<<root_ones_) - 1)) { }

    V root() const {
        return root_;
    }

    bool is_embeddable(const big_polymer<antichains_inst>& bp) const {
        return prefix_info(bp).dist_to_embeddable == 0;
    }

    int dist_to_embeddable(const big_polymer<antichains_inst>& bp) const {
        return prefix_info(bp).dist_to_embeddable;
    }

    GiNaC::ex rooted_embeddings(const big_polymer<antichains_inst>& bp) const {
        auto pi = prefix_info(bp);
        auto rl = mapping_delta_ + root_ones_;

        auto root_emb = GiNaC::binomial(n_, rl);
        if (n_ == 2*rl-2) {
            // Hack to make \binom{2d}{d-1} and \binom{2d}{d+1} print the same.
            // Clearly correct and may be removed if desired.
            root_emb = GiNaC::binomial(n_, n_ - rl);
        }

        return root_emb
            * GiNaC::binomial(rl, pi.changed_ones)
            * GiNaC::binomial(n_ - rl, pi.changed_zeros);
    }

    std::string debug_data(const big_polymer<antichains_inst>& bp) const {
        auto pi = prefix_info(bp);
        std::stringstream ss;
        ss << "root layer = " << mapping_delta_ + root_ones_
           << ", changed_zeros = " << pi.changed_zeros
           << ", changed_ones = " << std::to_string(pi.changed_ones)
           << ", rooted embeddings = " << rooted_embeddings(bp);
        return ss.str();
    }

    inline bool is_valid(const subpolymer<antichains_inst>& sp) const {
        assert(sp.mask_ != 0);

        int layer = sp.bp_.vtx(__builtin_ctz(sp.mask_)).layer();
        for (unsigned int it = sp.mask_; it; it &= it-1)
            if (layer != sp.bp_.vtx(__builtin_ctz(it)).layer())
                return false;

        return true;
    }

    /* Weight of polymer is is \lambda^|P| / (1+\lambda)^{|N(P) \cap M|}. */
    inline GiNaC::ex weight(const subpolymer<antichains_inst>& sp,
                            GiNaC::symbol lambda) const {
        auto pi = prefix_info(sp.bp_);
        assert(pi.dist_to_embeddable == 0);

        GiNaC::ex middle_layer_nghbrs = 0;
        std::set<V> nghbrs_prefix;

        for (int it = sp.mask_; it; it &= it-1) {
            V v = sp.bp_.vtx(__builtin_ctz(it));
            auto ambient_layer = mapping_delta_ + v.layer();

            int overcount = 0;
            for (unsigned int it2 = pi.special_coords; it2; it2 &= it2-1) {
                auto w = v.flip(__builtin_ctz(it2));
                if (w.layer() == middle_ones_) {
                    nghbrs_prefix.insert(w);
                    overcount++;
                }
            }

            if (v.layer() == middle_ones_ - 1)
                middle_layer_nghbrs += n_ - ambient_layer - overcount;
            else
                middle_layer_nghbrs += ambient_layer - overcount;
        }

        middle_layer_nghbrs += nghbrs_prefix.size();
        return GiNaC::pow(lambda, std::popcount(sp.mask_))
            / GiNaC::pow(1 + lambda, middle_layer_nghbrs);
    }
};

#endif
