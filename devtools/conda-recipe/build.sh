#!/bin/bash

# this step is done automatically by conda build
# cp -r $RECIPE_DIR/../.. $SRC_DIR
$PYTHON setup.py install

# Add more build steps here, if they are necessary.

# See
# http://docs.continuum.io/conda/build.html
# for a list of environment variables that are set during the build process.
