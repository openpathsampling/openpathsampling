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

Our coding style is mostly `PEP8`_ (with some allowances for things like
capitalization in variable names, when scientific custom prefers capitals),
and docstrings formatted according to the `numpydoc`_ standard.

.. _PEP8: https://www.python.org/dev/peps/pep-0008/
.. _numpydoc: https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt

-----

The sections linked below give details about extending OPS. The document on
code structure gives an overview of the code, and the other documents
describe specific aspects.

.. toctree::
    :maxdepth: 2

    code_structure
    pathsimulators
    engines
    pathmovers_and_movestrategies
    networks
    collective_variables
    storage_need_to_know
    
