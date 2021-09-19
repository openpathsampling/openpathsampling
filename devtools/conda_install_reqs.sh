#!/usr/bin/bash

# USAGE
# ./conda_install_reqs.sh
#   Install the requirements. If the environment variable OPS_ENV is set,
#   then create a new environment called $OPS_ENV. Note that you MUST set
#   the environment variable CONDA_PY (the version of Python to use, e.g.,
#   CONDA_PY="3.7").

DEVTOOLS_DIR=`dirname "${BASH_SOURCE[0]}"`

if [ ! command -v mamba ]
then
    EXE="conda"
else
    EXE="mamba"
fi

if [ ! -z "$OPS_ENV" ]
then
    $EXE create -q -y --name $OPS_ENV conda future pyyaml python=$CONDA_PY
    source activate $OPS_ENV
else
    $EXE install -y -q future pyyaml  # ensure that these are available
fi

# for some reason, these approaches to pinning don't always work (but conda
# always obeys if you explicitly request a pinned version)
#conda config --env --add pinned_packages python=$CONDA_PY
#cp $DEVTOOLS_DIR/../pinned $CONDA_PREFIX/

# WORKAROUNDS is normally empty; needed if other pkgs don't list all deps
WORKAROUNDS=""
REQUIREMENTS=`python ${DEVTOOLS_DIR}/setup_cfg_reqs.py`
TESTING=`python ${DEVTOOLS_DIR}/setup_cfg_reqs.py --extra test`
INTEGRATIONS=`cat ${DEVTOOLS_DIR}/tested_integrations.txt | tr "\n" " "`
EXPERIMENTAL=`cat ${DEVTOOLS_DIR}/experimental_reqs.txt | tr "\n" " "`
PY_INSTALL="python=$CONDA_PY"

INTEGRATIONS="future svgwrite ujson numpy scipy pandas matplotlib networkx netcdf4 psutil"

echo "PY_INSTALL=$PY_INSTALL"
echo "REQUIREMENTS=$REQUIREMENTS"
echo "INTEGRATIONS=$INTEGRATIONS"
echo "EXPERIMENTAL=$EXPERIMENTAL"
echo "WORKAROUNDS=$WORKAROUNDS"
echo "TESTING=$TESTING"

# TODO: allow different installs for different versions
# (needed this when msmbuilder was only available on older pythons; similar
# situations may come up in the future)
ALL_PACKAGES="$WORKAROUNDS " #$REQUIREMENTS" # $INTEGRATIONS" # $EXPERIMENTAL" # $TESTING"
ALL_PACKAGES="$WORKAROUNDS $EXPERIMENTAL " #$REQUIREMENTS" # $INTEGRATIONS" # $EXPERIMENTAL" # $TESTING"

echo "conda install -y -q -c conda-forge -c omnia $PY_INSTALL $ALL_PACKAGES"
$EXE install -y -q -c conda-forge -c omnia $PY_INSTALL $ALL_PACKAGES
