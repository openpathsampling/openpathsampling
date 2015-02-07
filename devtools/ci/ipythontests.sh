# Run ipython notebook tests

cd examples/ipython
testfail=0
python ipnbdoctest.py "alanine.ipynb" || testfail=1
python ipnbdoctest.py "sliced_sequential_ensembles.ipynb" || testfail=1
python ipnbdoctest.py "toy_dynamics_tis.ipynb" || testfail=1
python ipnbdoctest.py "toy_storage.ipynb" || testfail=1
python ipnbdoctest.py "langevin_integrator_check.ipynb" || testfail=1
python ipnbdoctest.py "sliced_sequential_ensembles.ipynb" || testfail=1
# python ipnbdoctest.py "visualization.ipynb" || testfail=1
cd ../..
if [ testfail -eq 1 ]
then
    exit 1
fi

