#!/bin/sh

echo travis_fold:start:nosetests
echo Run nose tests ...

testfail=0
pytest -vv --cov --cov-report xml:cov.xml || testfail=1
#nosetests -v -s --with-coverage || testfail=1
coveralls
echo travis_fold:end:nosetests

if [ $testfail -eq 1 ]
then
    exit 1
fi
