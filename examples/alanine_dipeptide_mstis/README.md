# Alanine dipeptide 6-state MSTIS

Multiple state transition interface sampling for alanine dipeptide. This
sequence is based on (although not identical to) the paper:

[W.-N. Du, K. A. Marino, and P. G. Bolhuis, “Multiple state transition interface sampling of alanine dipeptide in explicit solvent,” J. Chem. Phys., vol. 135, no. 14, p. 145102, 2011.](http://dx.doi.org/10.1063/1.3644344)

We obtain initial trajectories using the `FullBootstrapping` approach, and
then run replica exchange multiple state transition interface sampling.

- [`alanine_dipeptide_mstis_bootstrapping.ipynb`](http://github.com/openpathsampling/openpathsampling/blob/master/examples/ipython/alanine_dipeptide_mstis_bootstrapping.ipynb)
    [Part 1] Basic setup using 6 predefined states and their interfaces.
    This part will only construct initial trajectories that will be use
    later.
    
- [`alanine_dipeptide_mstis_prepare.ipynb`](http://github.com/openpathsampling/openpathsampling/blob/master/examples/ipython/alanine_dipeptide_mstis_bootstrapping.ipynb)
    [Part 2] This will use the previously generated paths and state
    definitions from a file and construct all necessary parts for the 6
    state MSTIS including the generation of initial minus paths for all
    states and equilibration. All results are stored in a file.

- [`alanine_dipeptide_mstis_run.ipynb`](http://github.com/openpathsampling/openpathsampling/blob/master/examples/ipython/alanine_dipeptide_mstis_run.ipynb)
    [Part 3] This will read the MSTIS setup from before and run a fixed
    steps and store the results. 

- [`alanine_dipeptide_mstis_restart.ipynb`](http://github.com/openpathsampling/openpathsampling/blob/master/examples/ipython/alanine_dipeptide_mstis_restart.ipynb)
    [Part 4] (Optional) This open the current production file and will run
    and append more steps.

- [`alanine_dipeptide_mstis_analysis.ipynb`](http://github.com/openpathsampling/openpathsampling/blob/master/examples/ipython/alanine_dipeptide_mstis_analysis.ipynb)
    [Part 5] Analysis of the results in the production file. Almost entirely a copy of the analysis of the toy MSTIS analysis notebook.
     


