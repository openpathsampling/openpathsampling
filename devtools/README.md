Developer Notes / Tools
=======================

Assorted notes for developers.

How to do a release
-------------------

1. Make a branch in the upstream `openpathsampling/openpathsampling`
   repository.
2. In that branch, change the version in `setup.cfg` to be a release version
   (not ending in `.dev0`, etc.)
3. Push that branch and make a PR against the `stable` branch. It should test,
   then perform a test deployment to `test.pypi`, then download the test
   deployment and run tests again. The logic for this is in the Autorelease
   repository. The text of the top post in this PR should be the release notes
   (use the write-release-notes script included in Autorelease to do this).
4. Merge the PR. This will create a GitHub release, and then take that release
   and push it to pypi. All of this logic is in Autorelease.


Docs Building & Hosting
-----------------------

After a travis build succeeds, the docs are built with sphinx and pushed to
the openpathsampling.org amazon s3 account (owned by John Chodera). The
credentials for that account are stored, encrypted, in the .travis.yml file.
(http://docs.travis-ci.com/user/build-configuration/#Secure-environment-variables)
