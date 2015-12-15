# This should run the nosetests in the right folder

# install python pip packages needed for testing etc.
# None of these must be necessary for build or installation
echo travis_fold:start:install.nose.pip.packages
echo Install pip packages
PIP_ARGS="-U"
$HOME/miniconda2/envs/${python}/bin/pip install $PIP_ARGS -r devtools/ci/nose-requirements-${python}.txt
echo travis_fold:end:install.nose.pip.packages

# install python conda packages needed for testing etc.
# None of these must be necessary for build or installation
echo travis_fold:start:install.nose.conda.packages
echo Install pip packages
$HOME/miniconda2/envs/${python}/bin/pip install $PIP_ARGS -r devtools/ci/nose-requirements-${python}.txt
echo travis_fold:end:install.nose.conda.packages

echo travis_fold:start:nosetests
echo Run nose tests ...

cd openpathsampling
cd tests
testfail=0
nosetests -v -s . || testfail=1
echo travis_fold:end:nosetests
cd ../..
if [ $testfail -eq 1 ]
then
    exit 1
fi
