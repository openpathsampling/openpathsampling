# Toy Model Multiple State Transition Interface Sampling

This example illustrates the use of OpenPathSampling for Multiple State
Transition Interface Sampling in a 2D toy potential with 3 states and using
a simple Langevin integrator.

- [`toy_mstis_1_setup.ipynb`](http://github.com/openpathsampling/openpathsampling/blob/master/examples/toy_model_mstis/toy_mstis_1_setup.ipynb)
    [Part 1] Setup the system and create initial pathways to be used in
    later parts. Contains an example on how to use bootstrapping.
- [`toy_mstis_2_run.ipynb`](http://github.com/openpathsampling/openpathsampling/blob/master/examples/toy_model_mstis/.ipynb)
    [Part 2] Use the previously generated initial pathways and generate
    data to be analyzed later.
- [`toy_mstis_3_analysis.ipynb`](http://github.com/openpathsampling/openpathsampling/blob/master/examples/toy_model_mstis/toy_mstis_3_analysis.ipynb)
    [Part 3] Take the previously generated data and do a complete analysis.
    Includes different visualizations, rate computations etc.
- [`toy_mstis_4_repex_analysis.ipynb`](http://github.com/openpathsampling/openpathsampling/blob/master/examples/toy_model_mstis/toy_mstis_4_repex_analysis.ipynb)
    [Part 4] Performs analysis of the replica exchange behavior, including
    analysis of replica mixing, replica flow, etc.
- [`toy_mstis_5_srtis.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/toy_model_mstis/toy_mstis_5_srtis.ipynb)
    [Part 5] Use the same network as used in `mstis.ipynb`, but use a
    single-replica TIS scheme. This notebook also includes some analysis of
    the SRTIS results.

## Appendix A: Splitting the analysis file

After running the calculations above, the following appendix files
demonstrate how to split the resulting file into two files, so that a file
without the full coordinate data can be used for analysis.

- [`toy_mstis_A1_split.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/toy_model_mstis/toy_mstis_A1_split.ipynb)
    [Part 1] Split the file into two.
- [`toy_mstis_A2_split_analysis.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/toy_model_mstis/toy_mstis_A2_split_analysis.ipynb)
    [Part 2] Perform analysis using the split file.
