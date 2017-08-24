#!/bin/bash
# This script was taken from https://github.com/pandegroup/mdtraj/tree/master/devtools

### Install Miniconda

echo travis_fold:start:install.conda
echo Install conda

if [ -z "$CONDA_PY" ]
then
    CONDA_PY=2.7
fi

pyV=${CONDA_PY:0:1}

MINICONDA=Miniconda${pyV}-latest-Linux-x86_64.sh
MINICONDA_MD5=$(curl -s https://repo.continuum.io/miniconda/ | grep -A3 $MINICONDA | sed -n '4p' | sed -n 's/ *<td>\(.*\)<\/td> */\1/p')
wget https://repo.continuum.io/miniconda/$MINICONDA
if [[ $MINICONDA_MD5 != $(md5sum $MINICONDA | cut -d ' ' -f 1) ]]; then
    echo "Miniconda MD5 mismatch"
    echo "Expected: $MINICONDA_MD5"
    echo "Found: $(md5sum $MINICONDA | cut -d ' ' -f 1)"
    exit 1
fi
bash $MINICONDA -b

export PATH=$HOME/miniconda${pyV}/bin:$PATH

# add omnia and update
conda config --add channels http://conda.anaconda.org/omnia
conda update --yes conda

echo travis_fold:end:install.conda
