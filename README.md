Assuming you have GiNaC installed, the following should create a `cluster` binary:

```
autoreconf -i
automake
./configure
make
```

TODO: Explain architectural decisions.

We also provide Jupyter notebooks for viewing coefficients (see below for specifics). Depending on your browser, you may not be able to horizontally scroll through the formulas. In this case, please open the Sage Jupyter notebook directly via Sage.

# Coefficients for `On Dedekind's problem, a sparse version of Sperner's theorem, and antichains of a given size in the Boolean lattice`

See [the Jupyter notebook](antichains.ipynb) to see a human-readable version of the coefficients, as well as details on how to consume the output of the program.

The `generated-antichains` directory contains pre-computed coefficients from the C++ program.

# Coefficients for `Counting independent sets in expanding bipartite regular graphs`

See [the Jupyter notebook](cartprod.ipynb) to see a human-readable version of the coefficients, as well as details on how to consume the output of the program.

The `generated-cartprod` directory contains pre-computed coefficients from the C++ program. The output file for the 6th coefficient was originally 44MiB in size, and for this reason it was simplified with Maxima's `ratsimp` function before being included in the repository.
