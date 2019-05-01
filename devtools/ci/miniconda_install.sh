#!/bin/bash
# This script was taken from https://github.com/pandegroup/mdtraj/tree/master/devtools

### Install Miniconda

if [ -z "$CONDA_PY" ]
then
    CONDA_PY=2.7
fi

pyV=${CONDA_PY:0:1}
conda_version="latest"
#conda_version="4.4.10"  # can pin a miniconda version like this, if needed

MINICONDA=Miniconda${pyV}-${conda_version}-Linux-x86_64.sh
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

# this puts the channel priority to (1) conda-forge; (2) omnia (3) defaults
conda config --add channels http://conda.anaconda.org/omnia
conda config --add channels http://conda.anaconda.org/conda-forge

# this may speed up conda's package resolution
conda config --set channel_priority strict

conda update --yes conda
