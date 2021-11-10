============
Installation
============

This page has instructions on installing OPS, installing software packages
that integrate with OPS, and testing that your OPS installation is working.

OpenPathSampling currently only works on POSIX systems (macOS and
UNIX/Linux). It is tested against Python 2.7, 3.7, 3.8, and 3.9.

Installing OpenPathSampling
===========================

.. _install-with-conda:

Installation with ``conda`` (recommended)
-----------------------------------------

OpenPathSampling can be installed with ``conda``, ``pip``, or by using
``setuptools``. We recommend using ``conda`` because it provides easy
installation of other tools, including molecular dynamics engines such as
OpenMM and GROMACS, which most users will want to install when using OPS.
Conda can be installed either with the `full Anaconda distribution
<https://www.anaconda.com/products/individual>`_, or with
the `smaller-footprint miniconda
<https://docs.conda.io/en/latest/miniconda.html>`_. 

OpenPathSampling can be installed from the ``conda-forge`` channel. To
install the most recent release, use the command ::

  $ conda install -c conda-forge openpathsampling

With that, you should be ready to use OPS!

Installation with ``pip`` or ``setuptools``
-------------------------------------------

OPS can also be installed with ``pip`` via the command ``python -m pip install
openpathsampling``, or the source can be downloaded and installed using
``setuptools`` and the command ``python setup.py install`` (run from the
root directory of the source repository). We only recommend these
installation mechanisms for users who are familiar with these tools.

Developer install
-----------------

OpenPathSampling can be installed as an "editable" (developer) install
through standard mechanisms (``python -m pip install -e .`` or ``python setup.py
develop``). 

Additionally, we provide a script that downloads the OPS source repository,
uses ``conda`` to install OPS's requirements (including several optional
integrations) and the installs OPS. Details on using this script are in
:ref:`quick-dev-install`.


Integration with other packages
===============================

When you install OpenPathSampling, you only get the core of OPS and the
internal toy dynamics engine. In practice, you'll probably want a more
powerful MD engine like OpenMM or GROMACS.

To add support for optional packages, all you need to do is install them.
For example, simply by having OpenMM or GROMACS installed on your system,
you will automatically have support for those engines within OPS. 

The tables below list several packages that provide OPS with extra
functionality when installed, and include links to their ``conda`` packages.

Molecular dynamics engines
--------------------------

+------------------------+-------------------------------------------------+
| | **OpenMM**           |                                                 |
| | |openmm-conda|_      |                                                 |
+------------------------+-------------------------------------------------+
| | **OpenMMTools**      | OPS has special support for integrators from    |
| | |openmmtools-conda|_ | OpenMMTools (and we strongly recommend using    |
|                        | one of the reversible integrators provided      |
|                        | there, such as VVVR).                           |
+------------------------+-------------------------------------------------+
| | **GROMACS**          | OPS will work with any installed GROMACS; it    |
| | |gromacs-conda|_     | does not need to come from ``conda``.           |
+------------------------+-------------------------------------------------+

.. |openmm-conda| image:: https://img.shields.io/conda/vn/conda-forge/openmm
.. |openmmtools-conda| image:: https://img.shields.io/conda/vn/conda-forge/openmmtools
.. |gromacs-conda| image:: https://img.shields.io/conda/vn/bioconda/gromacs

.. _openmm-conda: https://anaconda.org/conda-forge/openmm
.. _openmmtools-conda: https://anaconda.org/conda-forge/openmmtools
.. _gromacs-conda: https://anaconda.org/bioconda/gromacs

Collective variables
--------------------

In addition to the optional packages listed below, OPS ships with MDTraj
built-in.

+------------------------+-------------------------------------------------+
| | **PyEMMA**           |                                                 |
| | |pyemma-conda|_      |                                                 |
+------------------------+-------------------------------------------------+
| | **py-PLUMED**        |                                                 |
| | |plumed-conda|_      |                                                 |
+------------------------+-------------------------------------------------+

.. |pyemma-conda| image:: https://img.shields.io/conda/vn/conda-forge/pyemma
.. |plumed-conda| image:: https://img.shields.io/conda/vn/conda-forge/py-plumed

.. _plumed-conda: https://anaconda.org/conda-forge/py-plumed
.. _pyemma-conda: https://anaconda.org/conda-forge/pyemma


.. _run-tests:

Testing your installation
=========================

OpenPathSampling includes a thorough test suite, and running the test suite
is a good start to troubleshooting any installation problems. The OPS test
suite requires the packages ``pytest`` and (for legacy reasons) ``nose``.
These can be  installed with either ``conda`` or ``pip``. For example: ::

  $ conda install pytest nose

Once those are installed, you can run the test suite on your installation
with the command: ::

  $ py.test --pyargs openpathsampling.tests

The test suite includes over 1100 individual tests, and runs in around 2-3
minutes. All tests should either pass or skip.
