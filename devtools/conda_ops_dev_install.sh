#!/bin/sh

# This assumes you've already installed conda. Then it installs OPS as a
# developer install.

conda install -y -c conda-forge -c omnia openpathsampling
conda update -y --all
conda uninstall -y openpathsampling
git clone https://github.com/openpathsampling/openpathsampling.git
cd openpathsampling && pip install -e . && cd ..

# this puts to the OPS dev install in a directory so you can change to that
# directory and add remotes/switch branches as necessary
