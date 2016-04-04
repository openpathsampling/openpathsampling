.. _getting-started:

************
Installation
************

Platforms
=========

Mac and Linux. Windows only has limited OpenMM support. (Not recommended)


Install with Conda
------------------
.. _install-with-conda:

We recommend using `conda` to install OpenPathSampling.  `conda
<http://www.continuum.io/blog/conda>`_ is a python package manager built for
scientific python. Unlike ``easy_install`` or ``pip``, it handles binaries
and binary dependencies, which are critical for most scientific workflows.
If you're a ``conda`` user, you can install OpenPathSampling by adding the
omnia channel.

To install the most recent release of OpenPathSampling with conda, use the
following commands ::

  $ conda config --add channels omnia
  $ conda install openpathsampling

If you want the cutting edge of what OpenPathSampling can do, you can
install the development version using `conda install openpathsampling-dev`.
(Also, if you're reading this prior to the first official release.)

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

You will have to manually install the requirements and dependencies.

Dependencies
============

To use OpenPathSampling, the following libraries and software will need to
be installed. (TODO: this list is partial)

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
        Used for the examples and tutorials. Is optional, but highly
        recommended.


Testing Your Installation
=========================
Running the tests is a great way to verify that everything is working. The test
suite uses `nose <https://nose.readthedocs.org/en/latest/>`_, which you can pick
up via ``conda`` if you don't already have it. ::

  conda install nose

From the source directory ``openpathsampling/tests``, you can run the tests
by typing ``nosetests`` on the command line.
