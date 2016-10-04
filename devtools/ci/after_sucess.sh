#!/usr/bin/env bash
if [[ "$TRAVIS_PULL_REQUEST" != "false" ]]; then
    echo "This is a pull request. No deployment will be done."; exit 0
fi

if [[ "$TRAVIS_BRANCH" != "master" ]]; then
    echo "No deployment on BRANCH='$TRAVIS_BRANCH'"; exit 0
fi

PACKAGE_NAME=openpathsampling-dev
GROUP_NAME=omnia
BUILD_PATH=~/miniconda2/conda-bld

echo travis_fold:start:binstar.upload
if [[ "2.7" =~ "$python" ]]; then
    conda install --yes anaconda-client jinja2
    conda convert -p all ${BUILD_PATH}/linux-64/${PACKAGE_NAME}*.tar.bz2 -o ${BUILD_PATH}/
    anaconda -t ${ANACONDA_TOKEN}  upload  --force -u ${GROUP_NAME} -p ${PACKAGE_NAME} ${BUILD_PATH}/*/${PACKAGE_NAME}*.tar.bz2
fi

if [[ "$python" != "2.7" ]]; then
    echo "No deploy on PYTHON_VERSION=${python}"; exit 0
fi

echo travis_fold:end:binstar.upload

echo travis_fold:start:build.docs
# Create the docs and push them to S3
# -----------------------------------
conda install --yes pip
conda config --add channels http://conda.binstar.org/omnia
conda install --yes `conda build devtools/conda-recipe --output`
pip install numpydoc s3cmd
conda install --yes `cat docs/requirements.txt | xargs`

conda list -e

# Install pandoc for markdown support

(cd docs && python premake.py && make html && cd -)
ls -lt docs/_build
pwd
echo travis_fold:end:build.docs

echo travis_fold:start:upload.docs
python devtools/ci/push-docs-to-s3.py
echo travis_fold:end:upload.docs
