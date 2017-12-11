#!/bin/sh
# Run ipython notebook tests

ORIG_DIR=`pwd`
cd examples/
if [ "$CONDA_PY" = "2.7" ]; then
    ./ipynbtests.sh || exit 1
elif [ "$CONDA_PY" = "3.6" ]; then
    ./ipynbtests3.sh || exit 1
else
    echo "Notebook tests skipped for CONDA_PY=$CONDA_PY"
fi
cd $ORIG_DIR
