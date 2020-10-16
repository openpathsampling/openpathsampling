.. _AD-tps:

Flexible Length TPS on Alanine Dipeptide
========================================

This example shows how to set up, run, and analyze a 2-state TPS simulation.
The system studied is alanine dipeptide, which Peter Bolhuis likes to say is
"the hydrogen atom of biomolecules" -- it has enough complexity to
illustrate the ideas of biomolecular simulation, but is small enough that
even when solvated, a full TPS simulation can run on a laptop and be
converged overnight.

This example consists of three notebooks: one to create an initial
trajectory from a high-temperature run, one to run TPS, and one to analyze
the TPS results.  These notebooks can be found in the `OpenPathSampling
GitHub repository
<https://github.com/openpathsampling/openpathsampling/tree/master/examples/alanine_dipeptide_tps>`_.

.. toctree::

    alanine_dipeptide_tps/AD_tps_1_trajectory
    alanine_dipeptide_tps/AD_tps_2a_run_flex
    alanine_dipeptide_tps/AD_tps_3a_analysis_flex
