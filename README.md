Assuming you have GiNaC installed, the following should create a `cluster` binary:

```
autoreconf -i
automake
./configure
make
```

See [the Jupyter notebook](antichains-postprocessing.ipynb) to see a human-readable version of the coefficients, as well as details on how to consume the output of the program.

TODO: Explain architectural decisions.

# Generated outputs

The linked notebook above should be OK for viewing, but it contains big expressions.
To view full formulas, look at the exported sources in the `notebook-export` directory.
The `generated-antichains` directory contains pre-computed coefficients from the C++ program.
