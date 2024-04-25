#!/usr/bin/bash

# USAGE
# ./conda_install_reqs.sh
#   Install the requirements. If the environment variable OPS_ENV is set,
#   then create a new environment called $OPS_ENV. Note that you MUST set
#   the environment variable CONDA_PY (the version of Python to use, e.g.,
#   CONDA_PY="3.7").

DEVTOOLS_DIR=`dirname "${BASH_SOURCE[0]}"`

# use mamba if we can, otherwise use conda
if ! command -v mamba &> /dev/null; then
    EXE="conda"
else
    EXE="mamba"
fi

# If OPS_ENV is defined, we're creating a new environment. Otherwise install
# into current env.
if [ ! -z "$OPS_ENV" ]; then
    INSTALL="create --name $OPS_ENV"
else
    INSTALL="install"
fi

INSTALL_CMD="$EXE $INSTALL -y -q -c conda-forge -c omnia --override-channels"

# for some reason, these approaches to pinning don't always work (but conda
# always obeys if you explicitly request a pinned version)
#conda config --env --add pinned_packages python=$CONDA_PY
#cp $DEVTOOLS_DIR/../pinned $CONDA_PREFIX/

# WORKAROUNDS is normally empty; needed if other pkgs don't list all deps
WORKAROUNDS=""
REQUIREMENTS=`python ${DEVTOOLS_DIR}/setup_cfg_reqs.py`
TESTING=`python ${DEVTOOLS_DIR}/setup_cfg_reqs.py --extra test`
INTEGRATIONS=`cat ${DEVTOOLS_DIR}/tested_integrations.txt | tr "\n" " "`
EXPERIMENTAL=$(python ${DEVTOOLS_DIR}/setup_cfg_reqs.py --extra simstore)
PY_INSTALL="python=$CONDA_PY"

# PIP_INSTALLS is used for debugging installation problems -- override the
# default $REQUIREMENTS, etc. and move some installs to $PIP_INSTALLS
PIP_INSTALLS=""

echo "PY_INSTALL=$PY_INSTALL"
echo "REQUIREMENTS=$REQUIREMENTS"
echo "INTEGRATIONS=$INTEGRATIONS"
echo "EXPERIMENTAL=$EXPERIMENTAL"
echo "WORKAROUNDS=$WORKAROUNDS"
echo "TESTING=$TESTING"

# TODO: allow different installs for different versions
# (needed this when msmbuilder was only available on older pythons; similar
# situations may come up in the future)
ALL_PACKAGES="$WORKAROUNDS $REQUIREMENTS $INTEGRATIONS $EXPERIMENTAL $TESTING"

echo "$INSTALL_CMD $PY_INSTALL $ALL_PACKAGES"
if [ -z "$DRY" ]; then
    $INSTALL_CMD $PY_INSTALL $ALL_PACKAGES
fi

if [ -n "$OPS_ENV" ] && [ -z "$DRY" ]; then
    conda activate $OPS_ENV
fi

# occasional workaround; usually a do-nothing
if [ -n "$PIP_INSTALLS" ]; then
    echo "python -m pip install $PIP_INSTALLS"
    if [ -z "$DRY" ]; then
        python -m pip install $PIP_INSTALLS
    fi
fi
