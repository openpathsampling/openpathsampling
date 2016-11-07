# Toy Model Multiple Interface Set Transition Interface Sampling

This example illustrated the use of OpenPathSampling for [Multiple Interface
Set Transition Interface Sampling](http://dx.doi.org/10.1063/1.4890037) in a
2D toy model with 3 states and using a simple Langevin integrator.

- [`toy_mistis_1_setup_run.ipynb`](http://github.com/openpathsampling/openpathsampling/blob/master/examples/toy_model_mistis/toy_mistis_1_setup_run.ipynb)
    [Part 1] Setup the system and run MISTIS simulation.
- [`toy_mistis_2_flux.ipynb`](http://github.com/openpathsampling/openpathsampling/blob/master/examples/toy_model_mistis/toy_mistis_2_flux.ipynb)
    [Part 2; Optional] Calculate the flux for this system. Running this is
    optional because the result is provided in the next notebook; however,
    it would be needed to calculate the rate.
- [`toy_mistis_3_analysis.ipynb`](http://github.com/openpathsampling/openpathsampling/blob/master/examples/toy_model_mistis/toy_mistis_3_analysis.ipynb)
    [Part 3] Analyze the previously generated data to obtain, e.g., the rate
    matrix.
