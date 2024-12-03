#ifndef CLUSTER_CARTPROD_H
#define CLUSTER_CARTPROD_H

#include <iostream>
#include <vector>

template<typename B> class cartprod_vertex;
template<typename B> std::ostream &operator<<(std::ostream &,
                                              const cartprod_vertex<B> &);

template<typename B> class cartprod_vertex {
private:
    std::vector<typename B::V> data_;

public:
    cartprod_vertex(const std::vector<typename B::V>& coords)
        : data_(coords) { }

    unsigned int num_coords() const {
        return data_.size();
    }

    B::V get(unsigned int i) const {
        return data_[i];
    }

    cartprod_vertex<B> set(unsigned int i, B::V val) const {
        cartprod_vertex ret = *this;
        assert(ret.data_[i] != val);
        ret.data_[i] = val;
        return ret;
    }

    bool operator<(const cartprod_vertex<B>& cv) const {
        assert(cv.data_.size() == data_.size());
        return data_ < cv.data_;
    }

    unsigned int dist(const cartprod_vertex<B>& cv) const {
        unsigned int ans = 0;

        assert(cv.data_.size() == data_.size());
        for (int i = 0; i < data_.size(); i++)
            ans += data_[i].dist(cv.data_[i]);
        return ans;
    }

    unsigned int active_coords(const cartprod_vertex<B>& root) const {
        assert(data_.size() == root.data_.size());
        assert(data_.size() < 32);

        unsigned int ans = 0;
        for (int i = 0; i < data_.size(); i++)
            if (data_[i] != root.data_[i])
                ans |= 1<<i;

        return ans;
    }

    friend std::ostream& operator<< <B>(std::ostream& o,
                                        const cartprod_vertex<B>& bp);
};

template<typename B>
std::ostream& operator<<(std::ostream& o, const cartprod_vertex<B>& v) {
    o << "(" << v.data_[0];
    for (int i = 1; i < v.data_.size(); i++)
        o << ", " << v.data_[i];
    o << ")";
    return o;
}

/*
  A cartesian product of j copies of the same graph of size j.
*/
template<typename B> class cartprod {
private:
    B base_graph_;
    int j_;

public:
    using V = cartprod_vertex<B>;

    cartprod(int j) : base_graph_(B(j)), j_(j) { }

    std::vector<V> second_neighbourhood(V v) {
        assert(v.num_coords() == j_);
        std::vector<V> ans;

        for (int i = 0; i < j_; i++) {
            for (auto wi : base_graph_.neighbourhood(v.get(i))) {
                auto vi = v.set(i, wi);
                for (int j = i+1; j < j_; j++)
                    for (auto wj : base_graph_.neighbourhood(v.get(j))) {
                        auto vij = vi.set(j, wj);
                        ans.push_back(vij);
                    }
                }

            for (auto wii : base_graph_.second_neighbourhood(v.get(i))) {
                auto vii = v.set(i, wii);
                ans.push_back(vii);
            }
        }

        for (auto v2 : ans)
            assert(v.dist(v2) == 2);
        return ans;
    }
};

#endif
