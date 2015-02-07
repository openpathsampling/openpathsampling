# This should run the nosetests in the right folder

cd openpathsampling
cd tests
nosetests -v . || exit 1
cd ../..
