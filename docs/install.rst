.. _install:

************
Installation
************

OpenPathSampling currently only works on Mac and Linux. 

Install with Conda
==================
.. _install-with-conda:

We recommend using ``conda`` to install OpenPathSampling.  `conda
<http://www.continuum.io/blog/conda>`_ is a python package manager built for
scientific python. Unlike ``easy_install`` or ``pip``, it handles binaries
and binary dependencies, which are critical for most scientific workflows.
If you're a ``conda`` user, you can install OpenPathSampling by adding the
omnia channel.

To install the most recent release of OpenPathSampling with conda, use the
following commands ::

  $ conda config --add channels omnia
  $ conda install openpathsampling

.. note:: ``conda`` will automatically install all of the tricky dependencies
    from binary packages. This includes pytables, numpy, scipy, etc.  The
    easiest way to get conda is with the `Anaconda python distribution
    <https://store.continuum.io/cshop/anaconda/>`_.


Install from Source
===================

Clone the source code repository from github::

  $ git clone git://github.com/openpathsampling/openpathsampling.git

Then, in the directory containing the source code, you can install it with ::

  $ python setup.py install

You will have to manually install the requirements and dependencies.

Testing Your Installation
=========================
Running the tests is a great way to verify that everything is working. The test
suite uses `nose <https://nose.readthedocs.org/en/latest/>`_, which you can pick
up via ``conda`` if you don't already have it. ::

  conda install nose

From the source directory ``openpathsampling/tests``, you can run the tests
by typing ``nosetests`` on the command line.
