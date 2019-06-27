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

echo "Installing from ${REMOTE}/openpathsampling.git@$BRANCH"

git clone https://github.com/${REMOTE}/openpathsampling.git
pushd openpathsampling
git fetch
git checkout $BRANCH
git pull
source devtools/conda_install_reqs.sh
pip install -e .
popd
