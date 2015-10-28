## Python Notebook Examples
This folder contains some example python notebooks to learn and understand the ways OpenPathSampling works.

For this to work you need to have `ipython` and `ipython notebooks` installed. For help consult the 
[Installation Guide](http://ipython.org/install.html).

> Note there has been a major change in ipython notebooks which now run using jupyter server. 

Then open a terminal, walk to this folder (in your local copy) and run 
```
ipython notebook
```
or since version 4.0 you need to run a jupyter server using
```
jupyter notebook
```

## Examples
- Visualization [Notebook](http://nbviewer.ipython.org/github/jhprinz/msm-tis/blob/visualization/examples/ipython/Visualization%20Example.ipynb)
- Visualization [Notebook](http://nbviewer.ipython.org/github/jhprinz/msm-tis/blob/visualization/examples/ipython/Visualization%20Example.ipynb)

## Tutorials

Weina Alanine Example.ipynb [Notebook](http://nbviewer.ipython.org/github/jhprinz/msm-tis/blob/visualization/examples/Weina%20Alanine%20 Example.ipynb)

alanine.ipynb [Notebook](http://nbviewer.ipython.org/github/jhprinz/msm-tis/blob/visualization/examples/alanine.ipynb)

fast_sample_loading.ipynb [Notebook](http://nbviewer.ipython.org/github/jhprinz/msm-tis/blob/visualization/examples/fast_sample_loading.ipynb)

ipynbtest.py [Code](http://nbviewer.ipython.org/github/jhprinz/msm-tis/blob/visualization/examples/ipynbtest.py)

ipynbtest_tutorial.ipynb [Notebook](http://nbviewer.ipython.org/github/jhprinz/msm-tis/blob/visualization/examples/ipynbtest_tutorial.ipynb)
    An example notebook that explains how the ipython notebook testing script is used and what its features are. This is completely independent of _OpenPathSampling_ and it to be moved into a separate package.

- move_strategies_and_schemes.ipynb [Notebook](http://nbviewer.ipython.org/github/jhprinz/msm-tis/blob/visualization/examples/move_strategies_and_schemes.ipynb)


## Tests
- langevin_integrator_check.ipynb [Notebook](http://nbviewer.ipython.org/github/jhprinz/msm-tis/blob/visualization/examples/langevin_integrator_check.ipynb)
    A test notebook that checks that the toy_engine langevin integrator will actually sample from the correct distributions.

- mistis_setup.ipynb [Notebook](http://nbviewer.ipython.org/github/jhprinz/msm-tis/blob/visualization/examples/mistis_setup.ipynb)	
    [Part 1] of the MISTIS (Multi Interface Set TIS) testing notebooks. This creates the full setup and runs a few Monte Carlo steps for later analysis

- mistis_analysis.ipynb [Notebook](http://nbviewer.ipython.org/github/jhprinz/msm-tis/blob/visualization/examples/mistis_analysis.ipynb)	
    [Part 2] of the MISTIS (Multi Interface Set TIS) testing notebooks. This load the previously generated data, does some analysis and visualization on the results 



- mstis.ipynb [Notebook](http://nbviewer.ipython.org/github/jhprinz/msm-tis/blob/visualization/examples/mstis.ipynb)
    [Part 2] of the MSTIS (Multi State TIS) testing notebooks. This uses the previously generated initial pathways and generates data to be analyzed later

mstis_analysis.ipynb	Cleared output from notebooks	3 days ago
mstis_bootstrap.ipynb	Problem was in bootstrapping	a day ago
multistate_system_setup.ipynb	No wonder it was failing. Too many RETIS steps!	20 days ago
openmm_tutorial.ipynb	add openmm tutorial and remove folder	8 months ago
repex_networks.ipynb	Why doesn't Travis like this?	a day ago
sequential_ensembles.ipynb	fix nb	5 months ago
sliced_sequential_ensembles.ipynb	CVRangeVolume	5 months ago
storage_tutorial.ipynb	CVRangeVolume	5 months ago
toy_analysis.ipynb	Merge remote-tracking branch 'upstream/master' into engine_dt	5 months ago
toy_bootstrapping.ipynb	Descriptions on mstis_bootstrap	14 days ago
toy_plot_helpers.py	Live vis basically works!	4 months ago
toy_tis.ipynb	Using run_until	19 days ago
tree.pdf	deleted too large file	a year ago
troubleshooting_ops.ipynb	Added some draft documentation ipynbs	5 months ago
tutorial_visualization.ipynb	Merge remote-tracking branch 'upstream/master' into movescheme_minus_…	4 months ago
which_network.ipynb