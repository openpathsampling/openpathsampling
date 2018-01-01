#!/bin/sh
# Run ipython notebook tests

# TODO: completely remove this in favor of putting everything into
# examples/ipynbtests.sh

pushd examples/ && ./ipynbtests.sh || exit 1 && popd
