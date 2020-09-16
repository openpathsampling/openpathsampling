#!/usr/bin/bash

# USAGE
# ./conda_install_reqs.sh
#   Install the requirements. If the environment variable OPS_ENV is set,
#   then create a new environment called $OPS_ENV. Note that you MUST set
#   the environment variable CONDA_PY (the version of Python to use, e.g.,
#   CONDA_PY="3.7").

DEVTOOLS_DIR=`dirname "${BASH_SOURCE[0]}"`

#if [ ! -z "$OPS_ENV" ]
#then
    #conda create -q -y --name $OPS_ENV conda future pyyaml python=$CONDA_PY
    #source activate $OPS_ENV
#else
    #conda install -y -q future pyyaml  # ensure that these are available
#fi

# for some reason, these approaches to pinning don't always work (but conda
# always obeys if you explicitly request a pinned version)
#conda config --env --add pinned_packages python=$CONDA_PY
#cp $DEVTOOLS_DIR/../pinned $CONDA_PREFIX/

# WORKAROUNDS is normally empty; needed if other pkgs don't list all deps
WORKAROUNDS=""
REQUIREMENTS=`python ${DEVTOOLS_DIR}/setup_cfg_reqs.py`
TESTING=`python ${DEVTOOLS_DIR}/setup_cfg_reqs.py --extra testing`
INTEGRATIONS=`cat ${DEVTOOLS_DIR}/tested_integrations.txt | tr "\n" " "`
EXTRA=`cat ${DEVTOOLS_DIR}/optional_packages.txt | tr "\n" " "`
PY_INSTALL="python=$CONDA_PY"

echo "REQUIREMENTS=$REQUIREMENTS"
echo "WORKAROUNDS=$WORKAROUNDS"
echo "PY_INSTALL=$PY_INSTALL"
echo "TESTING=$TESTING"

# TODO: adjust this to be per python version
if [ "$CONDA_PY" != "3.7" ]; then
    ALL_PACKAGES="$WORKAROUNDS $REQUIREMENTS $TESTING $EXTRA"
else
    ALL_PACKAGES="$WORKAROUNDS $REQUIREMENTS $TESTING"  # no msmbuilder for py3.7?
fi

echo "conda install -y -q -c conda-forge -c omnia $PY_INSTALL $ALL_PACKAGES"
#conda install -y -q -c conda-forge -c omnia $PY_INSTALL $ALL_PACKAGES
