#!/usr/bin/bash

# REQUIRES ENV VARS: in particular, must define $CONDA_PY

# TODO: these should become cmd line arguments
DEVTOOLS_DIR=`dirname "${BASH_SOURCE[0]}"`

if [ ! -z "$OPS_ENV" ]
then
    conda create -q -y --name $OPS_ENV conda future pyyaml python=$CONDA_PY
    source activate $OPS_ENV
fi

# for some reason, these approaches to pinning don't always work (but conda
# always obeys if you explicitly request a pinned version)
#conda config --env --add pinned_packages python=$CONDA_PY
#cp $DEVTOOLS_DIR/../pinned $CONDA_PREFIX/
PACKAGES=`python ${DEVTOOLS_DIR}/install_recipe_requirements.py --dry ${DEVTOOLS_DIR}/conda-recipe/meta.yaml | tr "\n" " "`
TESTING=`cat ${DEVTOOLS_DIR}/testing_requirements.txt | tr "\n" " "`
EXTRA=`cat ${DEVTOOLS_DIR}/optional_packages.txt | tr "\n" " "`
PY_INSTALL="python=$CONDA_PY"
PINS=`cat ${DEVTOOLS_DIR}/../pinned | tr " " "" | tr "\n" " "`
ALL_PACKAGES="$PACKAGES $TESTING $EXTRA"
echo "conda install -y -q -c conda-forge -c omnia $PY_INSTALL $ALL_PACKAGES $PINS"
conda install -y -q -c conda-forge -c omnia $PY_INSTALL $ALL_PACKAGES $PINS
