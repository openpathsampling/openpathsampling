## Python Notebook Examples
This folder contains some example python notebooks to learn and understand the ways OpenPathSampling works.

For this to work you need to have `ipython` and `ipython notebooks` installed. For help consult the 
[Installation Guide](http://ipython.org/install.html).

Then open a terminal, walk to this folder `openpathsampling/examples/ipython` (in your local copy) and run
```
ipython notebook
```
or for version 4.0 and later
```
jupyter notebook
```
Note that starting with version 4.0 "the big split" the server running notebooks in kernels is separated from the actual ipython kernel, hence you run jupyter as the server for notebooks and need to have ipython installed in addition.

## Where should I start

There are a lot of different notebooks in here explaining (almost) every aspect of OpenPathSampling (OPS). There are a few typical notebooks to start with depending on your knowledge about PathSampling. 

If you are a _beginner_ and have not yet done any path sampling algorithms we recommend ...

If you are familiar with path sampling we still recommend starting with ...

See below for a complete list of available notebooks.

## Examples

#### MSTIS example (3 parts)

This is the main example illustrating the usage of OpenPathSampling for a Multi-State Transition Interface Sampling in a 2D toy potential with 3 states and using a simple Langevin integrator


- [`mstis_bootstrap.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/mstis_bootstrap.ipynb)
    [Part 1] of the MSTIS (Multi State TIS) testing notebooks. This will setup the general system and create initial pathways to be used in later parts. Contains an example on how to use bootstrapping.
- [`mstis.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/mstis.ipynb)
    [Part 2] of the MSTIS (Multi State TIS) testing notebooks. This uses the previously generated initial pathways and generates data to be analyzed later.
- [`mstis_analysis.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/mstis_analysis.ipynb)
    [Part 3] of the MSTIS (Multi State TIS) testing notebooks. This takes the previously generated data and does a complete analysis on them. Including different visualizations, rate computations, flow analysis, etc.

#### MISTIS example (2 parts)

This uses a similar setup as the MSTIS example but assignes different interfaces per outgoing transition. 

- [`mistis_setup.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/mistis_setup.ipynb)
    [Part 1] of the MISTIS (Multi Interface Set TIS) testing notebooks. This creates the full setup and runs a few Monte Carlo steps for later analysis

- [`mistis_analysis.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/mistis_analysis.ipynb)
    [Part 2] of the MISTIS (Multi Interface Set TIS) testing notebooks. This load the previously generated data, does some analysis and visualization on the results 

#### Alanine Dipeptide in explicit water example

- [`alanine.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/alanine.ipynb)
    A simple example noteboon on how to use OpenMM to run simulations on Alanine dipeptide in explicit solvent.

- [`Weina Alanine Example.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/Weina%20Alanine%20Example.ipynb)
    Basic setup from the W. Du and P. Bolhuis Paper on SRTIS. Is supposed to become the setup for an example in the publication.

## Tutorials

- [`fast_sample_loading.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/fast_sample_loading.ipynb)
    ... Seems to be removed in upcoming PRs

-  [`move_strategies_and_schemes.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/move_strategies_and_schemes.ipynb)
    A tutorial notebook to discuss the workings of move strategies and move schemes. Also a good starting point to understand how replica exchange moves are generated, detailed balance and all the little pitfalls there are about doing correct replica exchange moves.

- [`openmm_tutorial.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/openmm_tutorial.ipynb)
    Is a simple openmm tutorial originally written by A. Mey 

- [`repex_networks.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/repex_networks.ipynb)
    A test checking analysis functions for analysis of replica networks, i.e. treating the flow of replicas between ensembles as a graph and analyze it.

- [`sequential_ensembles.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/sequential_ensembles.ipynb)
    Incomplete tutorial on how to use the SequentialPathMover object.

- [`sliced_sequential_ensembles.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/sliced_sequential_ensembles.ipynb)
    ...

- [`storage_tutorial.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/storage_tutorial.ipynb)
    A tutorial notebook on how to work with the storage. Explains loading, saving, caching, etc...

- [`troubleshooting_ops.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/troubleshooting_ops.ipynb)
    This should become the FAQs of OpenPathSampling. A quick guide so solve the most common issues.

- [`tutorial_master.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/tutorial_visualization.ipynb)
    Contains a basic tutorial on visualization function. Still incomplete

- [`which_network.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/which_network.ipynb)
    A tutorial for beginners on how to decide which (MISTIS or MSTIS) to use. Still incomplete.

## Tests

- [`ipynbtest.py`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/ipynbtest.py)
    The python script to run and test ipython notebooks inside of travis. It is used to test the existing notebooks and use them as integration tests.

-  [`ipynbtest_tutorial.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/ipynbtest_tutorial.ipynb)
    An example notebook that explains how the ipython notebook testing script is used and what its features are. This is completely independent of _OpenPathSampling_ and it to be moved into a separate package.

- [`langevin_integrator_check.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/langevin_integrator_check.ipynb)
    A test notebook that checks that the toy_engine langevin integrator will actually sample from the correct distributions.

## Attic

- [`multistate_system_setup.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/multistate_system_setup.ipynb)
    Can be removed. Used to be part of the MSTIS example

- [`toy_bootstrapping.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/toy_bootstrapping.ipynb)
    [Part 1] of the large ToyDynamics example. Has been replaced by the mstis examples.

- [`toy_tis.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/toy_tis.ipynb)
    [Part 2] of the large ToyDynamics example. Has been replaced by the mstis examples.

- [`toy_analysis.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/toy_analysis.ipynb)
    [Part 3] of the large ToyDynamics example. Has been replaced by the mstis examples.


