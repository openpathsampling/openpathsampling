#!/usr/bin/bash

# USAGE
# ./conda_install_reqs.sh
#   Install the requirements. If the environment variable OPS_ENV is set,
#   then create a new environment called $OPS_ENV. Note that you MUST set
#   the environment variable CONDA_PY (the version of Python to use, e.g.,
#   CONDA_PY="3.7").

DEVTOOLS_DIR=`dirname "${BASH_SOURCE[0]}"`

if [ ! -z "$OPS_ENV" ]
then
    conda create -q -y --name $OPS_ENV conda future pyyaml python=$CONDA_PY
    source activate $OPS_ENV
else
    conda install -y -q future pyyaml  # ensure that these are available
fi

# for some reason, these approaches to pinning don't always work (but conda
# always obeys if you explicitly request a pinned version)
#conda config --env --add pinned_packages python=$CONDA_PY
#cp $DEVTOOLS_DIR/../pinned $CONDA_PREFIX/
PACKAGES=`python ${DEVTOOLS_DIR}/install_recipe_requirements.py --dry ${DEVTOOLS_DIR}/conda-recipe/meta.yaml | tr "\n" " "`
TESTING=`cat ${DEVTOOLS_DIR}/testing_requirements.txt | tr "\n" " "`
EXTRA=`cat ${DEVTOOLS_DIR}/optional_packages.txt | tr "\n" " "`
PY_INSTALL="python=$CONDA_PY"
PINS=`cat ${DEVTOOLS_DIR}/../pinned | tr -d " " | tr "\n" " "`
ALL_PACKAGES="$PACKAGES $TESTING $EXTRA"
echo "conda install -y -q -c conda-forge -c omnia $PY_INSTALL $ALL_PACKAGES $PINS"
conda install -y -q -c conda-forge -c omnia $PY_INSTALL $ALL_PACKAGES $PINS
