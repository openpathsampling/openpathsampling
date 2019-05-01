.. _install:

************
Installation
************

OpenPathSampling currently only works on Mac and Linux. For now, Python 2.7
is recommended, although most of the code also supports Python 3. Support
for Windows will be coming soon.

.. _install-with-conda:

Standard Install with Conda
===========================

We recommend using ``conda`` (with Python version 2.7) to install
OpenPathSampling.  `conda <http://www.continuum.io/downloads>`_ is a Python
package manager built for scientific Python, and which handles binary
dependencies seamlessly.  If you don't want the full ``conda`` installation,
the ``miniconda`` package provides much of the convenience of ``conda`` with
a smaller footprint.

OpenPathSampling is part of the ``omnia`` channel in ``conda``, although
some requirements are best found in the ``conda-forge`` channel. To install
the most recent release of OpenPathSampling with conda, use the following
command ::

  $ conda install -c conda-forge -c omnia openpathsampling

.. _developer-install-conda:

Developer Install with Conda
============================

To install a developer version of OPS (using ``conda``), change to a
directory where you want to OPS code (i.e., if you want the OPS git
repository at ``directory/openpathsampling/``, change to ``directory/``.)
Then download the ``conda_ops_dev_install.sh`` and run it ::

  $ curl -OLk https://raw.githubusercontent.com/openpathsampling/openpathsampling/master/devtools/conda_ops_dev_install.sh
  $ bash conda_ops_dev_install.sh

At this point, any changes to the code in that download of the OPS directory
will be live in your Python installation. You can use experimental code from
other forks by `adding the fork as a remote
<https://help.github.com/articles/adding-a-remote/>`_ and checking out a
branch.  You can combine changes from multiple users by merging them into a
branch in your local version of the repository.

.. _manual-install:

Manual Installation
===================

If you don't want to use ``conda``, you will have to manually obtain the
dependencies (advice on that coming soon). Then you can install from our
GitHub repository.

Clone the source code repository from github::

  $ git clone git://github.com/openpathsampling/openpathsampling.git

Then, in the directory containing the source code, you can install it with ::

  $ python setup.py install

Or, for a developer install ::

  $ python setup.py develop

Testing Your Installation
=========================
.. _run-tests:

Running the tests is a great way to verify that everything is working. The
test suite uses both `pytest <http://pytest.org>`_ and, for the time being,
`nose <https://nose.readthedocs.org/en/latest/>`_. You can pick these up via
``conda`` if you don't already have them. ::

  $ conda install pytest nose

From the source directory ``openpathsampling/tests``, you can run the tests
by typing ``py.test`` on the command line.

Beta testing experimental updates
=================================

In rare cases, you may want to test code that hasn't been merged into the
core of OPS yet. Instructions to install in this case are in the docs for
:ref:`install-devs`.
