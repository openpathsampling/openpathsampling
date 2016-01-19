#!/bin/sh

echo travis_fold:start:nosetests
echo Run nose tests ...

testfail=0
nosetests -v -s . || testfail=1
echo travis_fold:end:nosetests

if [ $testfail -eq 1 ]
then
    exit 1
fi
