# This should run the nosetests in the right folder

cd openpathsampling
cd tests
testfail=0
nosetests -v . || testfail=1
cd ../..
if [ testfail -eq 1 ]
then
    exit 1
fi
