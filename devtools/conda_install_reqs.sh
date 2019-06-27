#!/usr/bin/bash

# REQUIRES ENV VARS: in particular, must define $CONDA_PY

# TODO: these should become cmd line arguments
ENV_NAME="openpathsampling-py$CONDA_PY"
INCLUDE_PACKAGES="all"

DEVTOOLS_DIR=`dirname "${BASH_SOURCE[0]}"`

conda create -q -y --name $ENV_NAME conda future pyyaml python=$CONDA_PY
source activate $ENV_NAME

conda config --env --add pinned_packages python=$CONDA_PY
PACKAGES=`python ${DEVTO1OLS_DIR}/install_recipe_requirements.py --dry ${DEVTOOLS_DIR}/conda-recipe/meta.yaml | tr "\n" " "`
TESTING=`cat ${DEVTOOLS_DIR}/testing_requirements.txt | tr "\n" " "`
EXTRA=`cat ${DEVTOOLS_DIR}/optional_packages.txt | tr "\n" " "`
ALL_PACKAGES="$PACKAGES $TESTING $EXTRA"
echo $ALL_PACKAGES
conda install -y -q -c conda-forge -c omnia $ALL_PACKAGES
