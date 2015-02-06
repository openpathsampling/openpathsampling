# Run ipython notebook tests

cd examples/ipython
python ipnbdoctest.py "alanine.ipynb" || exit 1
python ipnbdoctest.py "sliced_sequential_ensembles.ipynb" || exit 1
python ipnbdoctest.py "toy_dynamics_tis.ipynb" || exit 1
python ipnbdoctest.py "toy_storage.ipynb" || exit 1
python ipnbdoctest.py "langevin_integrator_check.ipynb" || exit 1
python ipnbdoctest.py "sliced_sequential_ensembles.ipynb" || exit 1
# python ipnbdoctest.py "visualization.ipynb" || exit 1
cd ../..
