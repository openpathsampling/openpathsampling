#!/usr/bin/env bash
wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
bash Miniconda2-latest-Linux-x86_64.sh -b
export PATH=$HOME/miniconda2:$PATH
conda config --add channels omnia
conda install --yes openpathsampling

# execute script
