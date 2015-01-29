# Run ipython notebook tests

cd examples/ipython
python ipnbdoctest.py "alanine.ipynb"
python ipnbdoctest.py "sliced_sequential_ensembles.ipynb"
python ipnbdoctest.py "visualization.ipynb"
cd ../..