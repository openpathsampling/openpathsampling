echo $TRAVIS_PULL_REQUEST $TRAVIS_BRANCH

if [[ "$TRAVIS_PULL_REQUEST" != "false" ]]; then
    echo "This is a pull request. No deployment will be done."; exit 0
fi

if [[ "$TRAVIS_BRANCH" != "master" ]]; then
    echo "No deployment on BRANCH='$TRAVIS_BRANCH'"; exit 0
fi

if [[ "2.7" =~ "$python" ]]; then
    conda install --yes binstar jinja2
	conda convert -p all ~/miniconda/conda-bld/linux-64/openpathsampling-dev*.tar.bz2 -o ~/miniconda/conda-bld/
    binstar -t ${BINSTAR_TOKEN}  upload  --force --u omnia -p openpathsampling-dev $HOME/miniconda/conda-bld/*/openpathsampling-dev*.tar.bz2
fi

if [[ "$python" != "2.7" ]]; then
    echo "No deploy on PYTHON_VERSION=${python}"; exit 0
fi


# Create the docs and push them to S3
# -----------------------------------
conda install --yes pip
conda config --add channels http://conda.binstar.org/omnia
conda install --yes `conda build devtools/conda-recipe --output`
pip install numpydoc s3cmd
conda install --yes `cat docs/requirements.txt | xargs`

conda list -e

(cd docs && make html && cd -)
ls -lt docs/_build
pwd
python devtools/ci/push-docs-to-s3.py
