#!/usr/bin/env bash

# A script to test all the relevant IPython notebooks locally. We don't need
# to run this after every commit, but if you change much of the API, you
# should run it.

ipynbtest.py  alanine.ipynb
ipynbtest.py ipynbtest_tutorial.ipynb
ipynbtest.py langevin_integrator_check.ipynb
ipynbtest.py mistis_setup.ipynb
ipynbtest.py mistis_analysis.ipynb
ipynbtest.py move_strategies_and_schemes.ipynb
ipynbtest.py mstis_bootstrap.ipynb
ipynbtest.py mstis.ipynb
ipynbtest.py mstis_analysis.ipynb
ipynbtest.py repex_networks.ipynb
ipynbtest.py sequential_ensembles.ipynb
ipynbtest.py sliced_sequential_ensembles.ipynb
ipynbtest.py storage_mem_test.ipynb
ipynbtest.py storage_tutorial.ipynb
ipynbtest.py test_cv.ipynb
ipynbtest.py test_netcdfplus.ipynb
ipynbtest.py troubleshooting_ops.ipynb
ipynbtest.py tutorial_visualization.ipynb
ipynbtest.py which_network.ipynb
#ipynbtest.py Weina\ Alanine\ Example.ipynb
