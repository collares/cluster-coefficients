#include <iostream>
#include <ginac/ginac.h>
#include <string>
#include "antichains.h"
#include "hypercube.h"
#include "cartprod/cartprod_biclique_inst.h"

void print_coefficient(std::ostream& o,
                       GiNaC::ex expr,
                       std::string name,
                       GiNaC::symbol lambda) {
    o << name << ": " << expr;
    o << std::endl << std::endl;
    o << name << " at λ = 1: " << expr.subs(lambda == 1);
    o << std::endl << std::endl;
}

/* This is the main program output. It consists of two lines, each a function
   (but not a polynomial) in n and k. The first line corresponds to the number
   of rooted clusters in layer k-1 (it is multiplied by binomial(n, k-1)), and
   the second line corresponds to the number of rooted clusters in layer k+1
   (it is multiplied by binomial(n, k+1)).

   The output of the program will postprocessed in a Sage Jupyter notebook to
   make it a function of n, to evaluate it when λ = 1, and to compute other
   functions of the output which are used in the paper.
*/
void compute_antichains(int j, GiNaC::symbol lambda) {
    GiNaC::symbol n("n", "n"), k("k", "k");
    // When replacing n = 2*k, the expressions below should be equal.
    // This will be checked on the Jupyter notebook.
    auto Lj_lower = antichains_inst(-1, j, n, k).compute(lambda);
    std::cout << Lj_lower << std::endl;

    auto Lj_upper = antichains_inst(+1, j, n, k).compute(lambda);
    std::cout << Lj_upper << std::endl;

    auto Lj = Lj_lower + Lj_upper;
    print_coefficient(std::cerr, Lj, "L" + std::to_string(j), lambda);
}

/*
  The following function computes the coefficients for the "Counting
  independent sets in the hypercube revisited" paper.
*/
void compute_independent_sets(int j, GiNaC::symbol lambda) {
    GiNaC::symbol d("d", "d");
    auto Lj = hypercube_inst(j, d).compute(lambda);
    std::cout << Lj << std::endl;
    print_coefficient(std::cerr, Lj, "L" + std::to_string(j), lambda);
}

/*
  The following function computes the coefficients for the "Counting
  independent sets in bipartite regular cartesian product graphs" paper.
*/
void compute_cartprod_biclique(int j, GiNaC::symbol lambda) {
    GiNaC::symbol s("s", "s"), t("t", "t");
    auto Lj = cartprod_biclique_inst(j, s, t).compute(lambda);
    std::cout << Lj << std::endl;
    print_coefficient(std::cerr, Lj, "L" + std::to_string(j), lambda);
}

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "usage to compute L_j: " << argv[0] << " j" << std::endl;
        return 1;
    }

    int j = atoi(argv[1]);
    GiNaC::symbol lambda("λ", "\\lambda");
    compute_antichains(j, lambda);
}
