#!/usr/bin/bash

# REQUIRES ENV VARS: in particular, must define $CONDA_PY

# TODO: these should become cmd line arguments
ENV_NAME="openpathsampling-py$CONDA_PY"
DEVTOOLS_DIR=`dirname "${BASH_SOURCE[0]}"`

conda create -q -y --name $ENV_NAME conda future pyyaml python=$CONDA_PY
source activate $ENV_NAME

conda config --env --add pinned_packages python=$CONDA_PY
cp $DEVTOOLS_DIR/../pinned $CONDA_PREFIX/
PACKAGES=`python ${DEVTOOLS_DIR}/install_recipe_requirements.py --dry ${DEVTOOLS_DIR}/conda-recipe/meta.yaml | tr "\n" " "`
TESTING=`cat ${DEVTOOLS_DIR}/testing_requirements.txt | tr "\n" " "`
EXTRA=`cat ${DEVTOOLS_DIR}/optional_packages.txt | tr "\n" " "`
PY_INSTALL="python=$CONDA_PY"
#PINS=`cat ${DEVTOOLS_DIR}/../pinned | tr "\n" " "`
ALL_PACKAGES="$PACKAGES $TESTING $EXTRA"
echo "conda install -y -q -c conda-forge -c omnia $PY_INSTALL $ALL_PACKAGES" # $PINS"
conda install -y -q -c conda-forge -c omnia $PY_INSTALL $ALL_PACKAGES #$PINS
