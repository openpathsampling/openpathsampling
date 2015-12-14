# This should run the nosetests in the right folder

echo travis_fold:start:nosetests
echo Run nose tests ...

conda list -e
conda info
conda info -e

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
