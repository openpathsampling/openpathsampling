# Toy Model Multiple State Transition Interface Sampling

This example illustrates the use of OpenPathSampling for Multiple State
Transition Interface Sampling in a 2D toy potential with 3 states and using
a simple Langevin integrator.

- [`mstis_bootstrap.ipynb`](http://github.com/openpathsampling/openpathsampling/blob/master/examples/ipython/mstis_bootstrap.ipynb)
    [Part 1] Setup the system and create initial pathways to be used in
    later parts. Contains an example on how to use bootstrapping.
- [`mstis.ipynb`](http://github.com/openpathsampling/openpathsampling/blob/master/examples/ipython/mstis.ipynb)
    [Part 2] Use the previously generated initial pathways and generate
    data to be analyzed later.
- [`mstis_analysis.ipynb`](http://github.com/openpathsampling/openpathsampling/blob/master/examples/ipython/mstis_analysis.ipynb)
    [Part 3] Take the previously generated data and do a complete analysis.
    Includes different visualizations, rate computations, flow analysis,
    etc.
- [`srtis.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/srtis.ipynb)
    [Part 4] Use the same network as used in `mstis.ipynb`, but use a
    single-replica TIS scheme. This notebook also includes some analysis of
    the SRTIS results.
