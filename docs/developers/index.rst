.. _developer_info:

Overviews for Developers
========================

OpenPathSampling is an open-source project, and as such, we encourage
contributions from other developers.

Of course, OPS is also a rather large code, and getting into it can be
daunting. The documentation here is designed to help you identify what kinds
of objects you should be creating to make the improvements you want to make,
and to tell you what functions you must implement to make your additions
work in OPS.

If you're adding significant new functionality, we suggest a workflow where
you keep your additions in a new and separate repository, rather than in
your fork of OPS. This is because some functionality might not be applicable
to enough cases for us to include it in the core code. If your package is
outside OPS, it can still be distributed to those who are interested, even
if it isn't interoperable enough to go into the main repository.

.. rubric:: Coding Style

Our coding style is mostly `PEP8`_, with allowances for better readability
and matching standard scientific naming conventions (e.g., capitalization of
variable names). Docstrings are formatted according to the `numpydoc`_
standard.

.. _PEP8: https://www.python.org/dev/peps/pep-0008/
.. _numpydoc: https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt

.. rubric:: Tests and Examples

Additions to OPS should include thorough tests. We use `pytest`_ to run our
tests, although we require `nose`_ for some legacy portions of the tests.
Unit/integration tests are located in the ``openpathsampling/tests/``
directory.  Examples and integration tests should be written as Jupyter
notebooks and placed in the ``examples/`` directory.  Pull requests to the
OPS GitHub repository will automatically run all unit tests with each push.
We also run a subset of the system tests on each push, using `ipynbtest`_.

.. _pytest: http://pytest.org
.. _nose: http://nose.readthedocs.io/
.. _ipynbtest: https://github.com/jhprinz/ipynb-test

-----

The sections linked below give details about extending OPS. The document on
code structure gives an overview of the code, and the other documents
describe specific aspects.

.. toctree::
    :maxdepth: 2

    install_developers
    code_structure
    pathsimulators
    engines
    pathmovers_and_movestrategies
    networks
    collective_variables
    storage_need_to_know
    
