#!/usr/bin/bash

# REQUIRES ENV VARS: in particular, must define $CONDA_PY

# TODO: these should become cmd line arguments
ENV_NAME="openpathsampling-py$CONDA_PY"
DEVTOOLS_DIR=`dirname "${BASH_SOURCE[0]}"`

conda create -q -y --name $ENV_NAME conda future pyyaml python=$CONDA_PY
source activate $ENV_NAME

conda config --env --add pinned_packages python=$CONDA_PY
PACKAGES=`python ${DEVTOOLS_DIR}/install_recipe_requirements.py --dry ${DEVTOOLS_DIR}/conda-recipe/meta.yaml | tr "\n" " "`
TESTING=`cat ${DEVTOOLS_DIR}/testing_requirements.txt | tr "\n" " "`
EXTRA=`cat ${DEVTOOLS_DIR}/optional_packages.txt | tr "\n" " "`
PY_INSTALL="python=$CONDA_PY"
ALL_PACKAGES="$PACKAGES $TESTING $EXTRA"
#echo $ALL_PACKAGES
echo "conda install -y -q -c conda-forge -c omnia $PY_INSTALL $ALL_PACKAGES"
conda install -y -q -c conda-forge -c omnia $PY_INSTALL $ALL_PACKAGES
