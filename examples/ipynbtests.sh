#!/bin/sh
# Run ipython notebook tests

testfail=0

PYTHON_VERSION=`python -V 2>&1 | awk '{print $2}' | awk 'BEGIN { FS="." } { print $1 "." $2}'`
echo "Running tests for Python version: $PYTHON_VERSION"

dropbox_base_url="http://www.dropbox.com/s"

# TODO: slowly build this up as I get tests working for various test
# notebooks. Eventually, we get rid of the ipynbtests3 script and the
# version checking in devtools/ci/ipythontests.sh
case $PYTHON_VERSION in
    "2.7")
        mstis=$dropbox_base_url/1x4ny0c93gvu54n/toy_mstis_1k_OPS1.nc
        mistis=$dropbox_base_url/qaeczkugwxkrdfy/toy_mistis_1k_OPS1.nc
        ;;
    "3.5")
        mstis=$dropbox_base_url/1ulzssv5p4lr61f/toy_mstis_1k_OPS1_py36.nc
        mistis=$dropbox_base_url/8wldep8e26qignt/toy_mistis_1k_OPS1_py35.nc
        ;;
    "3.6")
        mstis=$dropbox_base_url/1ulzssv5p4lr61f/toy_mstis_1k_OPS1_py36.nc
        mistis=$dropbox_base_url/76981cbgxm639m3/toy_mistis_1k_OPS1_py36.nc
        ;;
    *)
        echo "Unsupported Python version: $PYTHON_VERSION"
esac

curl -OLk $mstis
curl -OLk $mistis
cp `basename $mstis` toy_mstis_1k_OPS1.nc
cp `basename $mistis` toy_mistis_1k_OPS1.nc

# from here should be the same for all versions
ls *nc
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
ipynbtest.py --strict --show-diff "test_openmm_integration.ipynb" || testfail=1
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
