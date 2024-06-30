#ifndef CLUSTER_HYPERCUBE_H
#define CLUSTER_HYPERCUBE_H

#include <bit>
#include <vector>
#include <ginac/ginac.h>
#include "cluster.h"
#include "hypercube_vertex.h"

class hypercube {
private:
    int sz_;

public:
    hypercube(int sz) : sz_(sz) { }

    std::vector<hypercube_vertex> second_neighbourhood(hypercube_vertex v) {
        std::vector<hypercube_vertex> ans;
        for (int i = 0; i < sz_; i++)
            for (int j = i+1; j < sz_; j++)
                ans.push_back(v.flip(i).flip(j));
        return ans;
    }
};

/*
  Morally this derives from instance, but we are using templates and
  monomorphisation instead. We could use the CRTP if desired, but I do not
  think we need it at the moment.
*/
class hypercube_inst : public instance<hypercube_inst> {
private:
    GiNaC::symbol d_;

public:
    using G = hypercube;
    using V = hypercube_vertex;
    hypercube graph_;

    // We change 2 coordinates at a time, so the hypercube must have at least
    // 2*j coordinates to capture every embeddable big polymer of size at most
    // j.
    hypercube_inst(int j, GiNaC::symbol d)
        : instance(j), graph_(hypercube(2*j)), d_(d) { }

    V root() const {
        return V(0);
    }

    unsigned int mask_union(const big_polymer<hypercube_inst>& bp) const {
        unsigned int ans = 0;
        for (int i = 0; i < bp.size(); i++)
            ans |= bp.vtx(i).active_coords(root());
        return ans;
    }

    std::string debug_data(const big_polymer<hypercube_inst>& bp) const {
        return "j = " + std::to_string(bp.size())
            + ", a = " + std::to_string(std::popcount(mask_union(bp)));
    }

    int dist_to_embeddable(const big_polymer<hypercube_inst>& bp) const {
        unsigned int mask = mask_union(bp);
        int d = std::popcount(~mask) - std::countl_zero(mask);
        return (d+1)/2;
    }

    int is_embeddable(const big_polymer<hypercube_inst>& bp) const {
        unsigned int mask = mask_union(bp);
        return !(mask & (mask+1));
    }

    GiNaC::ex rooted_embeddings(const big_polymer<hypercube_inst>& bp) const {
        int a = std::popcount(mask_union(bp));
        return GiNaC::pow(2, d_ - 1) * GiNaC::binomial(d_, a);
    }

    inline bool is_valid(const subpolymer<hypercube_inst>& sp) const {
        return true;
    }

    inline GiNaC::ex weight(const subpolymer<hypercube_inst>& sp,
                            GiNaC::symbol lambda) const {
        int a = std::popcount(mask_union(sp.bp_));

        GiNaC::ex num_nghbrs = std::popcount(sp.mask_);
        num_nghbrs *= d_ - a;

        std::set<V> nghbrs_prefix;
        for (int it = sp.mask_; it; it &= it-1) {
            int idx = __builtin_ctz(it);
            V v = sp.bp_.vtx(idx);
            for (int i = 0; i < a; i++)
                nghbrs_prefix.insert(v.flip(i));
        }

        num_nghbrs += nghbrs_prefix.size();
        return GiNaC::pow(lambda, std::popcount(sp.mask_))
            / GiNaC::pow(1 + lambda, num_nghbrs);
    }
};

#endif
