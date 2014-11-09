# This script was taken from https://github.com/pandegroup/mdtraj/tree/master/devtools

sudo apt-get update

### Install Miniconda

MINICONDA=Miniconda-latest-Linux-x86_64.sh
MINICONDA_MD5=$(curl -s http://repo.continuum.io/miniconda/ | grep -A3 $MINICONDA | sed -n '4p' | sed -n 's/ *<td>\(.*\)<\/td> */\1/p')
wget http://repo.continuum.io/miniconda/$MINICONDA
if [[ $MINICONDA_MD5 != $(md5sum $MINICONDA | cut -d ' ' -f 1) ]]; then
    echo "Miniconda MD5 mismatch"
    exit 1
fi
bash $MINICONDA -b

## Install conda pacakages

# This might make the --yes obsolete
# conda config --set always_yes yes --set changeps1 no

export PATH=$HOME/miniconda/bin:$PATH

conda update --yes conda
conda config --add channels http://conda.binstar.org/omnia

# Useful for debugging any issues with conda
conda info -a

conda create --yes -n ${python} --file devtools/ci/requirements-conda-${python}.txt
conda build devtools/conda-recipe
source activate $python
conda install --yes $HOME/miniconda/conda-bld/linux-64/yank-*

conda list -e

# install python pip packages

PIP_ARGS="-U"
$HOME/miniconda/envs/${python}/bin/pip install $PIP_ARGS -r devtools/ci/requirements-${python}.txt

# go back to the original directory we were in
# cd $MDTRAJ_DIR

pwd