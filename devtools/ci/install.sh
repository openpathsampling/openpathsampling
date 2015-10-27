# This script was taken from https://github.com/pandegroup/mdtraj/tree/master/devtools

sudo apt-get update

### Install Miniconda

echo travis_fold:start:install.conda
MINICONDA=Miniconda-latest-Linux-x86_64.sh
MINICONDA_MD5=$(curl -s https://repo.continuum.io/miniconda/ | grep -A3 $MINICONDA | sed -n '4p' | sed -n 's/ *<td>\(.*\)<\/td> */\1/p')
wget https://repo.continuum.io/miniconda/$MINICONDA
if [[ $MINICONDA_MD5 != $(md5sum $MINICONDA | cut -d ' ' -f 1) ]]; then
    echo "Miniconda MD5 mismatch"
    exit 1
fi
bash $MINICONDA -b
echo travis_fold:end:install.conda

## Install conda pacakages

# This might make the --yes obsolete
# conda config --set always_yes yes --set changeps1 no

export PATH=$HOME/miniconda/bin:$PATH

hash -r

echo travis_fold:start:install.conda.packages
conda config --add channels http://conda.anaconda.org/omnia
conda create --yes -n ${python} python=${python} --file devtools/ci/requirements-conda-${python}.txt
conda update --yes conda
source activate $python
echo travis_fold:end:install.conda.packages

# Useful for debugging any issues with conda
# conda info -a

echo travis_fold:start:install.pip.packages
# install python pip packages
PIP_ARGS="-U"
$HOME/miniconda/envs/${python}/bin/pip install $PIP_ARGS -r devtools/ci/requirements-${python}.txt
echo travis_fold:end:install.pip.packages
