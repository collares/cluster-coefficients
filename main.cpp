#include <iostream>
#include <ginac/ginac.h>
#include <string>
#include "antichains.h"
#include "hypercube.h"

void print_coefficient(std::ostream& o,
                       GiNaC::ex expr,
                       std::string name,
                       GiNaC::symbol lambda) {
    o << name << ": " << expr;
    o << std::endl << std::endl;
    o << name << " at λ = 1: " << expr.subs(lambda == 1);
    o << std::endl << std::endl;
}

void compute_antichains(int j,
                        GiNaC::symbol lambda,
                        GiNaC::symbol d) {
    for (int odd = 0; odd < 2; odd++) {
        GiNaC::ex n = 2*d + odd;
        std::string odd_str = odd ? "odd" : "even";
        std::cerr << "************* Case n = " << n << " *************"
                  << std::endl << std::endl;

        // Those are equal when n is even. We will test that on the Sage side.
        auto Lj_lower = antichains_inst(-1, j, n, d).compute(lambda);
        auto Lj_upper = antichains_inst(+1, j, n, d).compute(lambda);

        // To be piped into "sage antichains_postprocess.sage"
        std::cerr << "============ Start Sage input ============" << std::endl;
        std::cout << j << " " << odd << std::endl;
        std::cout << Lj_lower << std::endl;
        std::cout << Lj_upper << std::endl;
        std::cerr << "============  End  Sage input ============"
                  << std::endl << std::endl;

        auto Lj = Lj_lower + Lj_upper;
        print_coefficient(std::cerr,
                          Lj,
                          "L" + std::to_string(j) + " for " + odd_str + " n",
                          lambda);
    }
}

void compute_independent_sets(int j,
                              GiNaC::symbol lambda,
                              GiNaC::symbol d) {
    auto Lj = hypercube_inst(j, d).compute(lambda);
    print_coefficient(std::cerr, Lj, "L" + std::to_string(j), lambda);
}

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "usage to compute L_j: " << argv[0] << " j" << std::endl;
        return 1;
    }

    int j = atoi(argv[1]);
    GiNaC::symbol lambda("λ", "\\lambda"), d("d", "d");
    compute_antichains(j, lambda, d);
}
