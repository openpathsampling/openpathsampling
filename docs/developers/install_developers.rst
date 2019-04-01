.. _install-devs:

Installing for developers and beta-testers
==========================================

The overall process for adding code to OPS goes in 3 stages:

1. Code is developed in a user's fork of OPS.
2. Those changes are merged into the OPS master via a pull request on
   GitHub.
3. At some point, the master code is turned into a release.

For most users, code that has reached stage 3 is all that is needed. This is
what the standard install provides. Those who want the access to features
before they are released might be interested in code at stage 2. Installing
that is described in the "developer install" section of the installation
docs.

However, in rare cases, users might want to run the very raw code that is
only at stage 1. There are a few potential reasons for this. First, perhaps
you are developing new features, and want to use your own experimental code.
Second, perhaps there is an experimental feature that isn't ready for
general use, but you want to try it out very early. Third, perhaps you've
identified a bug, and we want to confirm with you that our fix does, in
fact, solve your problem. For any of these cases, you'll want to use the
installation procedure described here.

If you will be writing code for OPS, you should fork OPS and do the
developer install by cloning the repository of your own fork. For details,
see the `GitHub documentation on forking a repo
<https://help.github.com/en/articles/fork-a-repo>`_.

If you won't be contributing code (either you're beta-testing or confirming
a bugfix), you'll first do a :ref:`developer install
<developer-install-conda>`, then add a a "remote" so you can access code in
someone else's fork. You'll need to know which developer's fork you're
using, and what branch the code you'll be testing is on.

As an example, we'll add a remote from the developer |dwhswenson|_, and we
will call that remote ``dwhs`` (this name is up to you, but it is common to
use some short name such as the developer's initials). Then we'll check out a
branch called ``cool_new_feature`` from that fork, and give it the local
name ``dwhs_cool_new_feature``. Again, the local name is up to you, but
using the remote nickname and remote branch name make it easy to recognize.

.. |dwhswenson| replace:: ``dwhswenson``
.. _dwhswenson: http://github.com/dwhswenson

.. code::

    git remote add dwhs http://github.com/dwhswenson/openpathsampling.git
    git fetch dwhs
    git checkout -b dwhs_cool_new_feature dwhs/cool_new_feature

In the ``git checkout`` command, the first thing after the ``-b`` is the
local name of your branch, and the second is ``remote/branch_name``. 

With a developer installation (i.e., using ``pip install -e`` or ``setup.py
develop``), the code in the directory where you cloned the OPS repo is the
code you're running. So you can update to our more recent changes by running
``git pull``. You can switch to another version of the code (like a release
or the current master) by using ``git checkout`` to select the appropriate
branch or tag.

We strongly recommend that users working in this fashion become familiar
with git. Details like which branch you are in can make all the difference
when it comes to identifying problems.
