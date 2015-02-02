# Run ipython notebook tests

cd examples/ipython
python ipnbdoctest.py "alanine.ipynb"
python ipnbdoctest.py "sliced_sequential_ensembles.ipynb"
python ipnbdoctest.py "toy_dynamics_tis.ipynb"
python ipnbdoctest.py "toy_storage.ipynb"
python ipnbdoctest.py "langevin_integrator_check.ipynb"
python ipnbdoctest.py "sliced_sequential_ensembles.ipynb"
# python ipnbdoctest.py "visualization.ipynb"
cd ../..