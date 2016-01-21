#!/usr/bin/env sh

echo travis_fold:start:git_hash
pyversion=`python -c "import os; os.chdir(os.environ['HOME']); import openpathsampling as paths; print paths.version.full_version"`
gitversion=`git rev-parse HEAD`
echo travis_fold:end:git_has

echo "Installed version $pyversion from git commit hash $gitversion"
