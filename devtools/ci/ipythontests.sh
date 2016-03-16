#!/bin/sh
# Run ipython notebook tests

cd examples/ipython
testfail=0
#python ipynbtest.py "sliced_sequential_ensembles.ipynb" || testfail=1
date
ipynbtest.py "mstis_bootstrap.ipynb" || testfail=1
date
ipynbtest.py "mstis.ipynb" || testfail=1
date
ipynbtest.py "mstis_analysis.ipynb" || testfail=1
date
ipynbtest.py "repex_networks.ipynb" || testfail=1
date
ipynbtest.py "mistis_setup.ipynb" || testfail=1
date
ipynbtest.py "mistis_analysis.ipynb" || testfail=1
date
ipynbtest.py --strict "test_openmm_integration.ipynb" || testfail=1
date
ipynbtest.py "alanine.ipynb" || testfail=1

# needs to run after alanine since it need the trajectory.nc file
date
ipynbtest.py "tutorial_storage.ipynb" || testfail=1
date
ipynbtest.py --strict "test_netcdfplus.ipynb" || testfail=1
date
ipynbtest.py --strict "test_cv.ipynb" || testfail=1
date
ipynbtest.py --strict "test_pyemma.ipynb" || testfail=1

# python ipynbtest.py "visualization.ipynb" || testfail=1
cd ../..
if [ $testfail -eq 1 ]
then
    exit 1
fi
