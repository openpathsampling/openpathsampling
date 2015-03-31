.. _getting-started:

************
Installation
************

Platforms
=========

Mac and Linux. Windows only has limited OpenMM support. (Not recommended)

Supported Hardware
------------------

Needs OpenMM sofar!

Install with Conda
------------------
.. _install-with-conda:

`conda <http://www.continuum.io/blog/conda>`_ is a python package manager built for scientific python. Unlike ``easy_install`` or ``pip``, it handles binaries and binary dependencies, which are critical for most scientific workflows. If you're a ``conda`` user, you can install OpenPathSampling by adding the omnia channel. If you're not a conda user, you should look into it. ::

To install OpenPathSampling with conda, use the following commands ::

  $ conda config --add channels http://conda.binstar.org/omnia
  $ conda install openpathsampling

.. note:: ``conda`` will automatically install all of the tricky dependencies
    from binary packages automatically! This includes pytables / numpy / scipy!
    The easiest way to get conda is with the
    `Anaconda python distribution <https://store.continuum.io/cshop/anaconda/>`_.


Install from Source
-------------------
Clone the source code repository from github::

  $ git clone git://github.com/choderalab/openpathsampling.git

Then, in the directory containing the source code, you can install it with ::

  $ python setup.py install

Dependencies
============

To use OpenPathSampling, the following libraries and software will need to
be installed.

    Linux, Mac OS X or Windows operating system
        We develop mainly on 64-bit Linux and Mac machines. Windows is not
        well supported.

    `Python <http://python.org>`_ >= 2.7
        The development package (``python-dev`` or ``python-devel``
        on most Linux distributions) is recommended.

    `NumPy <http://numpy.scipy.org/>`_ >= 1.8.0
        Numpy is the base package for numerical computing in python.


Optional packages:

    `IPython <http://ipython.org>`_ >= 3.0.0
        Used for the examples and tutorials. Is optional, but highly recommended.


Testing Your Installation
=========================
Running the tests is a great way to verify that everything is working. The test
suite uses `nose <https://nose.readthedocs.org/en/latest/>`_, which you can pick
up via ``conda`` if you don't already have it. ::

  conda install nose

From the source directory ``openpathsampling/tests``, you can also run the tests
with ``nosetests .`` on the command line
