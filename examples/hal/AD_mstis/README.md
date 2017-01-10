# Alanine dipeptide 6-state MSTIS

Multiple state transition interface sampling for alanine dipeptide. This
sequence is based on (although not identical to) the paper:

[W.-N. Du, K. A. Marino, and P. G. Bolhuis, “Multiple state transition interface sampling of alanine dipeptide in explicit solvent,” J. Chem. Phys., vol. 135, no. 14, p. 145102, 2011.](http://dx.doi.org/10.1063/1.3644344)

We obtain initial trajectories using the `FullBootstrapping` approach, and
then run replica exchange multiple state transition interface sampling.

- [`AD_mstis_1_setup.ipynb`](http://github.com/openpathsampling/openpathsampling/blob/master/examples/alanine_dipeptide_mstis/AD_mstis_1_setup.ipynb)
    [Part 1] Basic setup using 6 predefined states and their interfaces.
    This part will only construct initial trajectories that will be use
    later.
    
- [`AD_mstis_2_run.ipynb`](http://github.com/openpathsampling/openpathsampling/blob/master/examples/alanine_dipeptide_mstis/AD_mstis_2_run.ipynb)
    [Part 2] This will read the MSTIS setup from before and run a fixed
    steps and store the results. 

- [`AD_mstis_3_restart.ipynb`](http://github.com/openpathsampling/openpathsampling/blob/master/examples/alanine_dipeptide_mstis/AD_mstis_3_restart.ipynb)
    [Part 3] (Optional) This open the current production file and will run
    and append more steps.

- [`AD_mstis_4_analysis.ipynb`](http://github.com/openpathsampling/openpathsampling/blob/master/examples/alanine_dipeptide_mstis/AD_mstis_4_analysis.ipynb)
    [Part 4] Analysis of the results in the production file. Almost entirely a copy of the analysis of the toy MSTIS analysis notebook.
     


