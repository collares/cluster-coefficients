#include <cassert>
#include <iostream>
#include <queue>
#include <set>
#include <stack>
#include <vector>
#include <ginac/ginac.h>

using namespace std;

const int NUM_PRINT_BITS = 6;

struct hypercube_vertex {
    unsigned int data;

    hypercube_vertex() : data(0) {
    }

    hypercube_vertex flip(int i) {
        hypercube_vertex n = *this;
        n.data ^= 1<<i;
        return n;
    }

    bool operator<(const hypercube_vertex& hv) const {
        return data < hv.data;
    }

    int dist(const hypercube_vertex& hv) {
        return __builtin_popcount(data ^ hv.data);
    }

    // Maybe implement set<hypercube_vertex> neighbourhood(int a);
};

ostream& operator<<(ostream& o, const hypercube_vertex& hv) {
    int v = hv.data;
    for (int i = 0; i < NUM_PRINT_BITS; i++, v /= 2)
        o << char('0' + (v%2));
    return o;
}

class big_polymer {
    // Annoyingly, during construction it is really convenient to
    // have this as a set but we want to use the ordering later.
    vector<hypercube_vertex> data;
    bool full_prefix;
    int prefix_sz;

public:
    big_polymer(const set<hypercube_vertex>& shv) {
        data = vector(shv.begin(), shv.end());

        int active_coords = 0;
        for (auto& v : data)
            active_coords |= v.data;

        active_coords += 1;
        full_prefix = !(active_coords & (active_coords-1));
        if (full_prefix)
            prefix_sz = __builtin_ctz(active_coords);

    }

    hypercube_vertex vtx(int idx) {
        return data[idx];
    }

    bool is_two_linked() const {
        assert(false);
    }

    bool operator<(const big_polymer& p) const {
        return data < p.data;
    }

    size_t num_vertices() const {
        return data.size();
    }

    pair<bool, int> active_coords_prefix() const {
        return {full_prefix, prefix_sz};
    }

    friend ostream& operator<<(ostream& o, const big_polymer& bp);
};

ostream& operator<<(ostream& o, const big_polymer& bp) {
    bool first = true;
    for (auto vtx : bp.data) {
        if (!first)
            o << ' ';
        o << vtx;
        first = false;
    }
    return o;
}

// Encodes a subset of an ordered big polymer.
struct subpolymer {
    big_polymer& bp;
    unsigned int mask;

    subpolymer(big_polymer& bp, unsigned int mask=0) : bp(bp), mask(mask) { }

    // Two polymers are incompatible if their union is 2-linked.
    bool is_incompatible_with(const subpolymer& other) const {
        for (int it1 = mask; it1; it1 &= it1-1) {
            hypercube_vertex v = bp.vtx(__builtin_ctz(it1));
            for (int it2 = other.mask; it2; it2 &= it2-1) {
                hypercube_vertex w = other.bp.vtx(__builtin_ctz(it2));
                if (v.dist(w) <= 2)
                    return true;
            }
        }

        return false;
    }

    size_t size() const {
        return __builtin_popcount(mask);
    }

    bool is_two_linked() const {
        assert(mask);

        queue<int> q;
        int visited; // also a mask
        int start = __builtin_ctz(mask);

        q.push(start);
        visited |= 1<<start;
        int comp_sz = 1;

        while (!q.empty()) {
            int idx = q.front();
            q.pop();

            for (int it = mask & ~visited; it; it &= it-1) {
                int nidx = __builtin_ctz(it);
                if (bp.vtx(idx).dist(bp.vtx(nidx)) == 2) {
                    q.push(nidx);
                    visited |= 1<<nidx;
                    comp_sz++;
                }
            }
        }

        return comp_sz == size();
    }

    GiNaC::ex weight(GiNaC::symbol lambda, GiNaC::symbol d) {
        int a = bp.active_coords_prefix().second;

        GiNaC::ex num_nghbrs = __builtin_popcount(mask);
        num_nghbrs *= (d-a);

        set<hypercube_vertex> nghbrs_prefix;
        for (int it = mask; it; it &= it-1) {
            int idx = __builtin_ctz(it);
            hypercube_vertex v = bp.vtx(idx);
            for (int i = 0; i < a; i++)
                nghbrs_prefix.insert(v.flip(i));
        }

        num_nghbrs += nghbrs_prefix.size();
        return pow(lambda, __builtin_popcount(mask))
            / pow(1 + lambda, num_nghbrs);
    }
};

ostream& operator<<(ostream& o, const subpolymer& sp) {
    bool first = true;
    o << "{";
    for (int it = sp.mask; it; it &= it-1) {
        int idx = __builtin_ctz(it);
        if (!first)
            o << ", ";
        o << sp.bp.vtx(idx);
        first = false;
    }
    return o << "}";
}

struct cluster {
    const big_polymer& bp;
    vector<subpolymer> elems;

    size_t total_sz;
    stack<unsigned int> union_mask; // to make this easily rewindable

    cluster(const big_polymer& bp) : bp(bp), total_sz(0) {
        assert(bp.active_coords_prefix().first);
        union_mask.push(0);
    }

    void push_subpolymer(const subpolymer& sp) {
        assert(&bp == &sp.bp);
        total_sz += sp.size();
        elems.push_back(sp);
        union_mask.push(union_mask.top() | sp.mask);
    }

    void pop_subpolymer() {
        total_sz -= elems.back().size();
        elems.pop_back();
        union_mask.pop();
    }

    bool covers_big_polymer() {
        return union_mask.top()+1 == 1<<bp.num_vertices();
    }

    size_t total_size() {
        return total_sz;
    }

    // The Ursell function of a graph is the sum over all subsets of
    // edges A inducing a connected graph of (-1)^{|A|}, normalized by
    // 1/v(H)!. This computes the Ursell function of the
    // incompatibility graph of the cluster.
    //
    // TODO: The Ursell function of a graph is, up to a multiplicative
    // factor of ((-1)^(1-v(H)))/v(H)!), an evaluation of the Tutte
    // polynomial at (1, 0). Tutte polynomials may be evaluated in
    // vertex-exponential time.
    //
    // TODO: The cluster is generated by adding polymers via backtracking.
    // There are many opportunities for reusing computation here.
    GiNaC::ex compute_ursell() {
        assert(elems.size());
        GiNaC::ex ans;

        vector<pair<int, int>> edges;
        for (int i = 0; i < elems.size(); i++)
            for (int j = i+1; j < elems.size(); j++)
                if (elems[i].is_incompatible_with(elems[j]))
                    edges.emplace_back(i, j);

        assert(edges.size() < 30);
        for (unsigned int mask = 0; mask < (1<<edges.size()); mask++) {
            vector<vector<int>> adj(elems.size());

            for (int it = mask; it; it &= it-1) {
                int i = __builtin_ctz(it);
                adj[edges[i].first].push_back(edges[i].second);
                adj[edges[i].second].push_back(edges[i].first);
            }

            queue<int> q;
            vector<bool> visited(elems.size());
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

            if (comp_sz != elems.size())
                continue;

            ans += __builtin_popcount(mask) % 2 ? -1 : 1;
        }

        return ans / GiNaC::factorial(elems.size());
    }

    GiNaC::ex compute_weight(GiNaC::symbol lambda, GiNaC::symbol d) {
        GiNaC::ex ans = compute_ursell();
        for (auto sp : elems)
            ans *= sp.weight(lambda, d);
        return ans;
    }

    friend ostream& operator<<(ostream& o, const cluster& c);
};

ostream& operator<<(ostream& o, const cluster& c) {
    bool first = true;
    o << "(";
    for (auto sp : c.elems) {
        if (!first)
            o << ", ";
        o << sp;
        first = false;
    }
    return o << ")";
}

// Generates L^a_m, the family of big polymers of size m whose active
// coordinates are precisely [a].
vector<big_polymer> generate_big_polymers(int active_prefix_sz /* a */,
                                          int num_elements /* m */,
                                          bool full_prefix) {
    assert(num_elements < 31);

    set<set<hypercube_vertex>> Lprev;
    set<set<hypercube_vertex>> Lcur;

    set<hypercube_vertex> start;
    start.insert(hypercube_vertex());
    Lprev.insert(start);

    /* Generate the lists L_m of 2-linked subsets of E of size m which
     * contain 0 and whose active coordinates are a subset of
     * [2*K]. */
    for (int k = 1; k < num_elements; k++) {
        Lcur.clear();
        for (const set<hypercube_vertex>& prev_set : Lprev)
            for (auto v : prev_set)
                for (int i = 0; i < active_prefix_sz; i++)
                    for (int j = 0; j < active_prefix_sz; j++)
                        if (i != j) {
                            set<hypercube_vertex> new_set = prev_set;
                            hypercube_vertex w = v.flip(i).flip(j);

                            if (new_set.find(w) == new_set.end()) {
                                new_set.insert(w);
                                assert(new_set.size() == k+1);
                                if (Lcur.find(new_set) == Lcur.end())
                                    Lcur.insert(new_set);
                            }
                        }

        swap(Lprev, Lcur);
    }

    vector<big_polymer> full_active;
    for (const set<hypercube_vertex>& sbp : Lprev) {
        big_polymer bp(sbp);
        if (!full_prefix ||
            bp.active_coords_prefix() == make_pair(true, active_prefix_sz))
            full_active.push_back(bp);
    }

    return full_active;
}

GiNaC::ex compute_Lk(GiNaC::symbol lambda, GiNaC::symbol d, int k) {
    GiNaC::ex ans;

    for (int m = 1; m <= k; m++) {
        GiNaC::ex a_sum;

        for (int a = 0; a <= 2 * k; a++) {
            GiNaC::ex acc;
            vector<big_polymer> Lam =
                generate_big_polymers(a, m, /* must_full_prefix = */ true);

            // Each set S in Lam ($L^a_m$) is potentially the vertex set of a
            // cluster. We will call such S a "big polymer". To obtain all
            // possible clusters, we will split each S into a tuple of
            // 2-linked subsets of total size k in all possible ways.

            // TODO: Convert to iterator
            for (auto &big_polymer : Lam) {
                cerr << "* Processing big polymer " << big_polymer
                     << " at j = " << m << ", a = " << a << "" << endl;

                vector<subpolymer> sps;

                assert(big_polymer.num_vertices() < 31);
                // Generate all possible 2-linked subpolymers
                // (polymers are always non-empty).
                //
                // TODO: There are cleverer ways to generate all
                // two-linked subsets. At the very least, one may keep
                // track of two-linked components.
                for (int msk = 1; msk < (1<<big_polymer.num_vertices()); msk++) {
                    subpolymer sp(big_polymer, msk);
                    if (sp.is_two_linked())
                        sps.push_back(sp);
                }

                // Now we want to find all ways of picking clusters,
                // which are tuples of subpolymers whose union is
                // big_polymer (possibly with repeats) and whose total
                // size is k.
                cluster current(big_polymer);

                function<void()> backtrack = [&]() {
                    if (current.total_size() > k)
                        return;

                    if (current.total_size() == k &&
                        current.covers_big_polymer()) {
                        GiNaC::ex cur = current.compute_weight(lambda, d);
                        cerr << "** Cluster " << current
                             << " has weight " << cur << endl;
                        acc += cur;
                    }

                    for (const subpolymer &sp : sps) {
                        current.push_subpolymer(sp);
                        backtrack();
                        current.pop_subpolymer();
                    }
                };

                backtrack();
            }

            a_sum += acc * GiNaC::binomial(d, a);
        }

        ans += a_sum/m;
    }

    return ans * pow(2, d-1);
}

int main() {
    int k;
    cin >> k;

    GiNaC::symbol lambda("L", "\\lambda"), d("d", "d");
    GiNaC::ex ans = compute_Lk(lambda, d, k);

    cout << "L_" << k << " = " << ans << endl;
    cout << "At lambda = 1: " << ans.subs(lambda == 1) << endl;
    cout << GiNaC::latex << ans << endl;
}
