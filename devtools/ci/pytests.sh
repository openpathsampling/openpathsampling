#!/bin/sh

echo travis_fold:start:pytests
echo Run pytest tests ...

testfail=0
pytest -vv -s --cov --cov-report xml:cov.xml || testfail=1
COVERALLS_PARALLEL=true coveralls
echo travis_fold:end:pytests

if [ $testfail -eq 1 ]
then
    exit 1
fi

if [ $CONDA_PY != "2.7" ] && [ -z "$MINIMAL" ]; then
    # experimental does not need to support Python 2
    echo travis_fold:start:experimental
    echo Running tests on experimental features
    pytest openpathsampling/experimental/ -vv -s || testfail=1
    echo travis_fold:end:experimental
fi


if [ $testfail -eq 1 ]
then
    exit 1
fi
