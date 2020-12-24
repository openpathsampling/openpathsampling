#/usr/bin/env python
from __future__ import print_function
import setup
import openpathsampling
from autorelease import DefaultCheckRunner, conda_recipe_version
from autorelease.version import get_setup_version
from packaging.version import Version

repo_path = '.'
SETUP_VERSION = get_setup_version(None, directory='.')
versions = {
    'package': openpathsampling.version.version,
    'netcdfplus': openpathsampling.netcdfplus.version.version,
    'setup.py': SETUP_VERSION,
}

RELEASE_BRANCHES = ['stable']
RELEASE_TAG = "v" + Version(SETUP_VERSION).base_version

if __name__ == "__main__":
    checker = DefaultCheckRunner(
        versions=versions,
        setup=setup,
        repo_path='.'
    )
    checker.release_branches = RELEASE_BRANCHES + [RELEASE_TAG]

    tests = checker.select_tests_from_sysargs()
    n_fails = checker.run_as_test(tests)
