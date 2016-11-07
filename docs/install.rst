.. _install:

************
Installation
************

OpenPathSampling currently only works on Mac and Linux, and requires Python
2.7. Support for Windows and Python 3 will be coming soon.

Install with Conda
==================
.. _install-with-conda:

We recommend using ``conda`` (with Python version 2.7) to install
OpenPathSampling.  `conda <http://www.continuum.io/downloads>`_ is a Python
package manager built for scientific Python, and which handles binary
dependencies seamlessly.  If you don't want the full ``conda`` installation,
the ``miniconda`` package provides much of the convenience of ``conda`` with
a smaller footprint.

OpenPathSampling is part of the ``omnia`` channel in ``conda``. To install
the most recent release of OpenPathSampling with conda, use the following
commands ::

  $ conda config --add channels omnia
  $ conda install openpathsampling


Install from Source
===================

If you would prefer to install from source, you will have to manually obtain
the dependencies (advice on that coming soon). Then you can install from our
GitHub repository.

Clone the source code repository from github::

  $ git clone git://github.com/openpathsampling/openpathsampling.git

Then, in the directory containing the source code, you can install it with ::

  $ python setup.py install

Testing Your Installation
=========================
Running the tests is a great way to verify that everything is working. The test
suite uses `nose <https://nose.readthedocs.org/en/latest/>`_, which you can pick
up via ``conda`` if you don't already have it. ::

  conda install nose

From the source directory ``openpathsampling/tests``, you can run the tests
by typing ``nosetests`` on the command line.
