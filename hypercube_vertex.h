#ifndef CLUSTER_HYPERCUBE_VERTEX_H
#define CLUSTER_HYPERCUBE_VERTEX_H

#include <bit>

namespace params {
    const int num_print_bits = 8;
}

struct hypercube_vertex {
    unsigned int data_;

    hypercube_vertex(unsigned int data) : data_(data) {}

    hypercube_vertex flip(int i) {
        hypercube_vertex n = *this;
        n.data_ ^= 1<<i;
        return n;
    }

    bool operator<(const hypercube_vertex& hv) const {
        return data_ < hv.data_;
    }

    int dist(const hypercube_vertex& hv) const {
        return std::popcount(data_ ^ hv.data_);
    }

    int layer() const {
        return std::popcount(data_);
    }

    unsigned int active_coords(const hypercube_vertex& root) const {
        return data_ ^ root.data_;
    }
};

std::ostream& operator<<(std::ostream& o, const hypercube_vertex& v) {
    int d = v.data_;
    for (int i = 0; i < params::num_print_bits; i++, d /= 2)
        o << char('0' + (d%2));
    return o;
}

#endif
