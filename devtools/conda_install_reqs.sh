#!/usr/bin/bash

# USAGE
# ./conda_install_reqs.sh
#   Install the requirements. If the environment variable OPS_ENV is set,
#   then create a new environment called $OPS_ENV. Note that you MUST set
#   the environment variable CONDA_PY (the version of Python to use, e.g.,
#   CONDA_PY="3.7").

DEVTOOLS_DIR=`dirname "${BASH_SOURCE[0]}"`

if ! command -v mamba &> /dev/null
then
    EXE="conda"
else
    EXE="mamba"
fi

INSTALL_CMD="$EXE install -y -q -c conda-forge -c omnia --override-channels"

# TODO: is this preinstall needed? We're certainly no longer using pyyaml,
# and it looks like we should be, but aren't, using future in setup_cfg_reqs
if [ ! -z "$OPS_ENV" ]
then
    $INSTALL_CMD --name $OPS_ENV conda future pyyaml python=$CONDA_PY
    source activate $OPS_ENV
else
    $INSTALL_CMD future pyyaml  # ensure that these are available
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


REQUIREMENTS="future svgwrite numpy scipy pandas networkx matplotlib psutil
ujson mdtraj netcdf4"
PIPREQS=""
INTEGRATIONS=""
EXPERIMENTAL=""

# PIP_INSTALLS is used for debugging installation problems -- override the
# default $REQUIREMENTS, etc. and move some installs to $PIP_INSTALLS
PIP_INSTALLS="${PIPREQS}"
#TESTING=""

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
$INSTALL_CMD $PY_INSTALL $ALL_PACKAGES

if [ -n "$PIP_INSTALLS" ]
then
    echo "python -m pip install $PIP_INSTALLS"
    python -m pip install $PIP_INSTALLS
fi
