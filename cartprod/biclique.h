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
  This is a complete bipartite graph of size sz_. The bipartition consists of
  evens and odds.
*/
class biclique {
private:
    int sz_;

public:
    using V = biclique_vertex;

    biclique(int sz) : sz_(sz) { }

    std::vector<V> neighbourhood(V v) {
        std::vector<V> ret;

        int cur = !v.parity();
        while (cur < sz_) {
            ret.push_back(cur);
            cur += 2;
        }
        return ret;
    }

    std::vector<V> second_neighbourhood(V v) {
        std::vector<V> ret;

        int cur = v.parity();
        while (cur < sz_) {
            if (cur != v.data_)
                ret.push_back(cur);
            cur += 2;
        }
        return ret;
    }
};

#endif
