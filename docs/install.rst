.. _install:

************
Installation
************

OpenPathSampling currently only works on Mac and Linux. It is tested against
Python 2.7, 3.6, and 3.7, although there may be some corners of the code
that aren't Python 3-compatible yet.

.. note:: As of OpenPathSampling 1.1, OpenMM will no longer be automatically
          installed when you install OPS. However, the OpenMM engine will be
          immediately available if you install OpenMM yourself. See the
          `OpenMM installation instructions
          <http://docs.openmm.org/latest/userguide/application.html#installing-openmm>`_
          for a detailed guide, but ``conda install -c conda-forge -c omnia
          openmm`` will work for most people. (COMING SOON: details on OPS
          integrations with other tools.)

.. _install-with-conda:

Standard Install with Conda
===========================

We recommend using ``conda`` to install OpenPathSampling.  `conda
<http://www.continuum.io/downloads>`_ is a Python package manager built for
scientific Python, and which handles binary dependencies seamlessly.  If you
don't want the full ``conda`` installation, the ``miniconda`` package
provides much of the convenience of ``conda`` with a smaller footprint.

OpenPathSampling is part of the ``conda-forge`` channel in ``conda``.  To
install the most recent release of OpenPathSampling with conda, use the
following command ::

  $ conda install -c conda-forge openpathsampling

With that, you should be ready to use OPS!

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

Additional functionality of the ``conda_ops_dev_install.sh`` script is
described in :ref:`quick-dev-install`.

.. _manual-install:

Manual Installation
===================

If you don't want to use ``conda``, you will have to manually obtain the
dependencies, which you can see listed under ``install_requires`` in
``setup.py``. Then you can install from our GitHub repository.

Clone the source code repository from GitHub::

  $ git clone https://github.com/openpathsampling/openpathsampling.git

Then, in the directory containing the source code, you can install it with ::

  $ python setup.py install

Or, for a developer install ::

  $ python setup.py develop

.. _run-tests:

Testing Your Installation
=========================

Running the tests is a great way to verify that everything is working. The
test suite uses `pytest <http://pytest.org>`_ and, for legacy reasons, also
requires the `nose <https://nose.readthedocs.org/en/latest/>`_ package. You can pick these up via ``conda`` if you don't already have them. ::

  $ conda install pytest nose

From the source directory ``openpathsampling/tests``, you can run the tests
by typing ``py.test`` on the command line. The test suite includes over 900
individual tests, and runs in around 2-3 minutes.

Beta testing experimental updates
=================================

In rare cases, you may want to test code that hasn't been merged into the
core of OPS yet. Instructions to install in this case are in the docs for
:ref:`install-devs`.
