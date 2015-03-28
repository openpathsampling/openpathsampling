#!/bin/sh
# Run ipython notebook tests

cd examples/ipython
testfail=0
python ipynbtest.py "alanine.ipynb" || testfail=1
python ipynbtest.py "sliced_sequential_ensembles.ipynb" || testfail=1
python ipynbtest.py "toy_boostrapping.ipynb" || testfail=1
python ipynbtest.py "toy_tis.ipynb" || testfail=1
python ipynbtest.py "toy_storage.ipynb" || testfail=1
python ipynbtest.py "langevin_integrator_check.ipynb" || testfail=1
python ipynbtest.py "sliced_sequential_ensembles.ipynb" || testfail=1
# needs to run after alanine since it need the trajectory.nc file
python ipynbtest.py "storage_tutorial.ipynb" || testfail=1
# python ipynbtest.py "visualization.ipynb" || testfail=1
cd ../..
if [ $testfail -eq 1 ]
then
    exit 1
fi

