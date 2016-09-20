#!/bin/sh
# Run ipython notebook tests

testfail=0
#python ipynbtest.py "sliced_sequential_ensembles.ipynb" || testfail=1
cd toy_model_mstis/
date
ipynbtest.py "toy_mstis_1_setup.ipynb" || testfail=1
date
ipynbtest.py "toy_mstis_2_run.ipynb" || testfail=1
date
# TODO: curl down toy_mstis_v1.nc
ipynbtest.py "toy_mstis_3_analysis.ipynb" || testfail=1
date
#ipynbtest.py "srtis.ipynb" || testfail=1
#date
ipynbtest.py "toy_mstis_4_repex_analysis.ipynb" || testfail=1
cd ../toy_model_mistis/
date
ipynbtest.py "toy_mistis_1_setup_run.ipynb" || testfail=1
date
# skip toy_mistis_2_flux: not needed
# TODO: curl down toy_mistis_v1.nc
ipynbtest.py "toy_mistis_3_analysis.ipynb" || testfail=1
date
cd ../tests/
cp ../toy_model_mstis/mstis.nc ./
ipynbtest.py --strict --show-diff "test_openmm_integration.ipynb" || testfail=1
date
ipynbtest.py --strict --show-diff "test_snapshot.ipynb" || testfail=1
date
ipynbtest.py --strict --show-diff "test_netcdfplus.ipynb" || testfail=1
date
ipynbtest.py --strict --show-diff "test_cv.ipynb" || testfail=1
date
ipynbtest.py --strict --show-diff "test_pyemma.ipynb" || testfail=1
date
cd ../misc/
cp ../toy_model_mstis/mstis.nc ./
ipynbtest.py "tutorial_storage.ipynb" || testfail=1

cd ../..
if [ $testfail -eq 1 ]
then
    exit 1
fi
