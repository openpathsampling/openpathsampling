# This script was taken from https://github.com/pandegroup/mdtraj/tree/master/devtools

sudo apt-get update

### Install Miniconda

echo travis_fold:start:install.conda
echo Install conda

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

export PATH=$HOME/miniconda2/bin:$PATH

hash -r

echo travis_fold:start:install.conda.packages
echo Install conda packages

conda config --add channels http://conda.anaconda.org/omnia
conda create --yes -n ${python} python=${python}
conda update --yes conda
source activate $python
echo travis_fold:end:install.conda.packages