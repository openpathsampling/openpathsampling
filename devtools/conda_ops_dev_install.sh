#!/bin/sh

# This assumes you've already installed conda. Then it installs OPS as a
# developer install.

# USAGE:
# ./conda_ops_dev_install.sh
#    install the master branch of the openpathsampling (main) fork
# ./conda_ops_dev_install.sh REMOTE
#    install the master branch of the fork REMOTE
# ./conda_ops_dev_install.sh REMOTE BRANCH
#    install BRANCH from the fork REMOTE

# defaults
REMOTE="openpathsampling"
BRANCH="master"

if [ $# -eq 2 ]; then
    REMOTE=$1
    BRANCH=$2
elif [ $# -eq 1 ]; then
    REMOTE=$1
elif [ $# -gt 2 ]; then
    echo "Invalid number of arguments: $# ($@)"
    exit 1
fi

if [ -z "$CONDA_PY" ]; then
    # if undefined, use major/minor of current Python
    CONDA_PY=`python -V 2>&1 |    # redirect stderr to stdout
              head -n 1 |         # ensure we only take first line
              cut -d " " -f2 |    # get the version number
              cut -d "." -f1,2`   # only keep major/minor
fi

if [ -z "$OPS_ENV" ]; then
    ENV_NAME="environment $CONDA_DEFAULT_ENV"
else
    ENV_NAME="new environment $OPS_ENV"
fi

echo "Installing from ${REMOTE}/openpathsampling.git@$BRANCH"
echo "Installing with Python $CONDA_PY in $ENV_NAME"

git clone https://github.com/${REMOTE}/openpathsampling.git
pushd openpathsampling
git fetch
git checkout $BRANCH
git pull
source devtools/conda_install_reqs.sh
pip install --no-deps -e .
popd
