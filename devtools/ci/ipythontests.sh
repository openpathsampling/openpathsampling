#!/bin/sh
# Run ipython notebook tests

ORIG_DIR=`pwd`
cd examples/
./ipynbtests.sh || exit 1
cd $ORIG_DIR
