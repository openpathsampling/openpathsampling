#!/bin/sh

echo travis_fold:start:pytests
echo Run pytest tests ...

testfail=0
pytest -vv -s --cov --cov-report xml:cov.xml || testfail=1
#nosetests -v -s --with-coverage || testfail=1
coveralls
echo travis_fold:end:pytests

if [ $testfail -eq 1 ]
then
    exit 1
fi
