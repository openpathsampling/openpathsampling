.. beta-testing:

Beta testing experimental updates
=================================

The overall process for adding code to OPS goes something like this:

1. Code is developed in a user's fork of OPS.
2. Those changes are merged into the OPS master.
3. At some point, the master code is turned into a release.

For most users, the release (which can be installed by conda), is all you
need. If you want to live on he bleeding edge, you can use the current
master, which is what you get when if you perform a developer install. You
can update this by running the command ``git pull``.

In rare cases, you might want to use code that hasn't been merged into the
OPS master yet. The primary cases for that are either (a) you plan to
develop for OPS; (b) you want to test a highly experimental feature; or (c)
we're asking you to check whether a bugfix does, in fact, fix a bug you've
reported.

If you will be writing code for OPS, you should fork OPS and do the
developer install by checking out the repository of your own fork. For
details, see the GitHub documentation.

In the other cases, you'll first do a developer install, then add what's
called a "remote" so you can access code in someone else's fork. You'll need
to know which user's fork you're using, and what branch the code you'll be
testing in on.

As an example, we'll add a remote from the developer |dwhswenson|_, and we
will call that remote ``dwhs`` (this name is up to you, but it is common to
use some short name such as the developers initials). Then we'll check out a
branch called ``cool_new_feature`` from that fork, and give it the local
name ``dwhs_cool_new_feature``. Again, the local name is up to you, but
using the remote nickname and remote branch name make it easy to recognize.

.. |dwhswenson| replace:: ``dwhswenson``
.. _dwhswenson: http://github.com/dwhswenson

.. code::

    git remote add dwhs http://github.com/dwhswenson/openpathsampling
    git fetch dwhs
    git checkout -b dwhs_cool_new_feature dwhs/cool_new_feature

In the ``git checkout`` comamnd, the first thing after the ``-b`` is the
local name of your branch, and the second is ``remote/branch_name``. 

