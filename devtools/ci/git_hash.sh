#!/usr/bin/env sh

# NOTE: the pyversion is determined by ensuring we're not in the git repo
# directory, where the local openpathsampling directory (without a
# version.py) is the version we load
pyversion=`python -c "import os; os.chdir(os.environ['HOME']); import openpathsampling as paths; print paths.version.full_version"`
gitversion=`git rev-parse HEAD`
echo "Installed version $pyversion from git commit hash $gitversion"
