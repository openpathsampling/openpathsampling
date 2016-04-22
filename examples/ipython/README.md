## Python Notebook Examples
This folder contains some example python notebooks to learn and understand
the ways OpenPathSampling works.

For this to work you need to have `ipython` and `ipython notebooks`
installed. For help consult the [Installation
Guide](http://ipython.org/install.html).

Then open a terminal, walk to this folder
`openpathsampling/examples/ipython` (in your local copy) and run
```
ipython notebook
```
or for version 4.0 and later
```
jupyter notebook
```
Note that starting with version 4.0 "the big split" the server running
notebooks in kernels is separated from the actual ipython kernel, hence you
run jupyter as the server for notebooks and need to have ipython installed
in addition.

## Where should I start

There are a lot of different notebooks in here explaining (almost) every
aspect of OpenPathSampling (OPS). There are a few typical notebooks to start
with depending on your knowledge about PathSampling. 

If you are a _beginner_ and have not yet done any path sampling algorithms
we recommend ...

If you are familiar with path sampling we still recommend starting with ...

See below for a complete list of available notebooks.

## Examples

#### Simple TPS examples

Transition path sampling (TPS) is the simplest path sampling algorithm.
These examples demostrate just how easy it is to set up a TPS simulation
using OpenPathSampling.

* [`simple_tps.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/simple_tps.ipynb)
    This example does flexible-pathlength TPS. You need to start with (a) a
    working MD engine for your system; (b) state definitions for your
    system; and (c) a trajectory (not necessarily physically real, but as
    close as possible) connecting the initial and final states for your
    system.
* [`simple_fixed_length_tps.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/simple_fixed_length_tps.ipynb)
    This example is for fixed-pathlength TPS. In addition to the information
    from the flexible-pathlength version, you also need to know an
    appropriate pathlength to use: here we use 5000 frames.

#### MSTIS example (3 parts)

This is the main example illustrating the usage of OpenPathSampling for a
Multi-State Transition Interface Sampling in a 2D toy potential with 3
states and using a simple Langevin integrator


- [`mstis_bootstrap.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/mstis_bootstrap.ipynb)
    [Part 1] of the MSTIS (Multi State TIS) testing notebooks. This will
    setup the general system and create initial pathways to be used in later
    parts. Contains an example on how to use bootstrapping.
- [`mstis.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/mstis.ipynb)
    [Part 2] of the MSTIS (Multi State TIS) testing notebooks. This uses the
    previously generated initial pathways and generates data to be analyzed
    later.
- [`mstis_analysis.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/mstis_analysis.ipynb)
    [Part 3] of the MSTIS (Multi State TIS) testing notebooks. This takes
    the previously generated data and does a complete analysis on them.
    Including different visualizations, rate computations, flow analysis,
    etc.

#### MISTIS example (2 parts)

This uses a similar setup as the MSTIS example but assignes different
interfaces per outgoing transition. 

- [`mistis_setup.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/mistis_setup.ipynb)
    [Part 1] of the MISTIS (Multi Interface Set TIS) testing notebooks. This
    creates the full setup and runs a few Monte Carlo steps for later
    analysis

- [`mistis_analysis.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/mistis_analysis.ipynb)
    [Part 2] of the MISTIS (Multi Interface Set TIS) testing notebooks. This
    load the previously generated data, does some analysis and visualization
    on the results 

#### Alanine Dipeptide in explicit water example

- [`alanine.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/alanine.ipynb)
    A simple example noteboon on how to use OpenMM to run simulations on
    Alanine dipeptide in explicit solvent.

- [`alanine_dipeptide_old_example.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/alanine_dipeptide_old_example.ipynb)
    Basic setup from the W. Du and P. Bolhuis Paper on SRTIS. Is supposed to
    become the setup for an example in the publication.

#### Alanine Weina Du 6 state MSTIS setting

- [`alanine_dipeptide_mstis_bootstrapping.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/alanine_dipeptide_mstis_bootstrapping.ipynb)
    [Part 1] Basic setup from the W. Du and P. Bolhuis (2011) Paper on SRTIS using 6 predefined states and their interfaces. This part will only construct initial trajectories that will be use later.[^1]
    
[^1]: W.-N. Du, K. A. Marino, and P. G. Bolhuis, “Multiple state transition interface sampling of alanine dipeptide in explicit solvent,” J. Chem. Phys., vol. 135, no. 14, p. 145102, 2011.
    
- [`alanine_dipeptide_mstis_prepare.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/alanine_dipeptide_mstis_bootstrapping.ipynb)
	[Part 2] This will use the previously generated paths and state definitions from a file and construct all necessary parts for the 6 state MSTIS including the generation of initial minus paths for all states and equilibration. All results are stored in a file.

- [`alanine_dipeptide_mstis_run.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/alanine_dipeptide_mstis_run.ipynb)
    [Part 3] This will read the MSTIS setup from before and run a fixed steps and store the results. 

- [`alanine_dipeptide_mstis_restart.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/alanine_dipeptide_mstis_restart.ipynb)
    [Part 4] This open the current production file and will run and append more steps.

- [`alanine_dipeptide_mstis_analysis.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/alanine_dipeptide_mstis_analysis.ipynb)
    [Part 5] Analysis of the results in the production file. Almost entirely a copy of the analysis of the toy MSTIS analysis notebook.
     

## Tutorials

-  [`move_strategies_and_schemes.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/move_strategies_and_schemes.ipynb)
    A tutorial notebook to discuss the workings of move strategies and move
    schemes. Also a good starting point to understand how replica exchange
    moves are generated, detailed balance and all the little pitfalls there
    are about doing correct replica exchange moves.

- [`repex_networks.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/repex_networks.ipynb)
    A test checking analysis functions for analysis of replica networks,
    i.e. treating the flow of replicas between ensembles as a graph and
    analyze it.

- [`committors.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/committors.ipynb)
    Examples of running committor calculations, as well as doing a committor
    analysis based on shooting points from a path sampling simulation that
    has already been run.

- [`sequential_ensembles.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/sequential_ensembles.ipynb)
    Incomplete tutorial on how to use the SequentialPathMover object.

- [`sliced_sequential_ensembles.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/sliced_sequential_ensembles.ipynb)
    ...

- [`tutorial_storage.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/tutorial_storage.ipynb)
    A tutorial notebook on how to work with the storage. Explains loading,
    saving, caching, etc...

- [`troubleshooting_ops.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/troubleshooting_ops.ipynb)
    This should become the FAQs of OpenPathSampling. A quick guide so solve
    the most common issues.

- [`tutorial_visualization.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/tutorial_visualization.ipynb)
    Contains a basic tutorial on visualization function. Still incomplete

- [`which_network.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/which_network.ipynb)
    A tutorial for beginners on how to decide which (MISTIS or MSTIS) to
    use. Still incomplete.

## Tests

Many of these notebooks are also used at integration tests. We use the
[`ipynbtest.py`](http://github.com/jhprinz/ipynb-test) script to run
notebooks from the command line, and, when appropriate, check the
correctness of the results.

- [`localtest.sh`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/localtest.sh)
    A test script to run *all* the notebooks listed here. Only a subset are
    used as continuous integration tests. Requires
    [`ipynbtest.py`](http://github.com/jhprinz/ipynb-test).

- [`langevin_integrator_check.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/langevin_integrator_check.ipynb)
    A test notebook that checks that the toy_engine langevin integrator will
    actually sample from the correct distributions.

- [`test_snapshot_modifier.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/test_snapshot_modifier.ipynb)
    A notebook to test and visually inspect whether the `RandomVelocities`
    snapshot modifier gives the correct distribution.

## Attic

- [`toy_bootstrapping.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/toy_bootstrapping.ipynb)
    [Part 1] of the large ToyDynamics example. Has been replaced by the mstis examples.

- [`toy_tis.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/toy_tis.ipynb)
    [Part 2] of the large ToyDynamics example. Has been replaced by the mstis examples.

- [`toy_analysis.ipynb`](http://github.com/choderalab/openpathsampling/blob/master/examples/ipython/toy_analysis.ipynb)
    [Part 3] of the large ToyDynamics example. Has been replaced by the mstis examples.
