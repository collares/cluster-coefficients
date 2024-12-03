#ifndef CLUSTER_BICLIQUE_H
#define CLUSTER_BICLIQUE_H

struct biclique_vertex {
    unsigned int data_;

    biclique_vertex(unsigned int data) : data_(data) {}

    bool operator<(const biclique_vertex& bv) const {
        return data_ < bv.data_;
    }

    bool operator!=(const biclique_vertex& bv) const {
        return data_ != bv.data_;
    }

    bool operator==(const biclique_vertex& bv) const {
        return data_ == bv.data_;
    }

    unsigned int parity() const {
        return data_ & 1;
    }

    // This shouldn't be here. dist should be a method of the graph.
    int dist(const biclique_vertex& bv) const {
        if (*this == bv)
            return 0;
        if (parity() != bv.parity())
            return 1;
        return 2;
    }
};

std::ostream& operator<<(std::ostream& o, const biclique_vertex& v) {
    o << v.data_;
    return o;
}

/*
  This is a K_{j, j}, where j = sz/2. The bipartition consists of evens and
  odds.
*/
class biclique {
private:
    int j_;

public:
    using V = biclique_vertex;

    biclique(int sz) : j_(sz/2) { }

    std::vector<V> neighbourhood(V v) {
        unsigned int parity_bit = !v.parity();

        std::vector<V> ret;
        for (int i = 0; i < j_; i++)
            ret.push_back((2*i) | parity_bit);
        return ret;
    }

    std::vector<V> second_neighbourhood(V v) {
        unsigned int parity_bit = v.parity();

        std::vector<V> ret;
        for (int i = 0; i < j_; i++)
            if (i != v.data_/2)
                ret.push_back((2*i) | parity_bit);
        return ret;
    }
};

#endif
