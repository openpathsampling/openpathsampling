#!/bin/sh
# Run ipython notebook tests

testfail=0
curl -OLk http://www.dropbox.com/s/1x4ny0c93gvu54n/toy_mstis_1k_OPS1.nc
curl -OLk http://www.dropbox.com/s/qaeczkugwxkrdfy/toy_mistis_1k_OPS1.nc
ls *nc
#python ipynbtest.py "sliced_sequential_ensembles.ipynb" || testfail=1
cd toy_model_mstis/
date
ipynbtest.py "toy_mstis_1_setup.ipynb" || testfail=1
date
ipynbtest.py "toy_mstis_2_run.ipynb" || testfail=1
date
ipynbtest.py "toy_mstis_3_analysis.ipynb" || testfail=1
date
ipynbtest.py "toy_mstis_4_repex_analysis.ipynb" || testfail=1
#date
#ipynbtest.py "toy_mstis_5_srtis.ipynb" || testfail=1
cd ../toy_model_mistis/
date
ipynbtest.py "toy_mistis_1_setup_run.ipynb" || testfail=1
date
# skip toy_mistis_2_flux: not needed
ipynbtest.py "toy_mistis_3_analysis.ipynb" || testfail=1
date
cd ../tests/
cp ../toy_model_mstis/mstis.nc ./
ipynbtest.py --strict "test_openmm_integration.ipynb" || testfail=1
date
ipynbtest.py --strict "test_snapshot.ipynb" || testfail=1
date
ipynbtest.py --strict "test_netcdfplus.ipynb" || testfail=1
date
ipynbtest.py --strict "test_cv.ipynb" || testfail=1
date
ipynbtest.py --strict "test_pyemma.ipynb" || testfail=1
date
cd ../misc/
cp ../toy_model_mstis/mstis.nc ./
ipynbtest.py "tutorial_storage.ipynb" || testfail=1

cd ..
rm toy_mstis_1k_OPS1.nc
rm toy_mistis_1k_OPS1.nc
if [ $testfail -eq 1 ]
then
    exit 1
fi
