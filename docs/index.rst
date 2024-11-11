OpenPathSampling
================

A Python library to facilitate path sampling algorithms. 

OpenPathSampling (OPS) makes it easy to perform many variants of transition
path sampling (TPS) and transition interface sampling (TIS), as well as
other useful calculations for rare events, such as committor analysis and
flux calculations. In addition, it is a powerful library to build new path
sampling methods.

OPS is independent of the underlying molecular dynamics engine, and
currently has support for OpenMM and Gromacs, as well as an internal engine
suitable for 2D toy models.

To learn more about what OPS can do, look at our :ref:`examples <examples>`.
If you want to jump right in, take a look at how easy it is to
:ref:`install <install>`!

OPS is an open-source project, distributed under the MIT license. Join us in
the development process on GitHub_, and follow `@pathsampling
<http://twitter.com/pathsampling>`_ on Twitter for updates!

.. _GitHub: http://github.com/openpathsampling/openpathsampling

.. note:: Documentation is still in progress. Please see :ref:`getting_help`
          for how to contact us with questions.

To see the most recent updates to the code, see the `release notes
<https://github.com/openpathsampling/openpathsampling/releases>`_ page on
GitHub.

Citing
------

OPS was described in a pair of papers published in JCTC:

1. David W.H. Swenson, Jan-Hendrik Prinz, Frank Noé, John D. Chodera, and
   Peter G. Bolhuis. "OpenPathSampling: A flexible, open framework for path
   sampling simulations. 1. Basics." J. Chem. Theory Comput. **15**, 813
   (2019).
   https://doi.org/10.1021/acs.jctc.8b00626
2. David W.H. Swenson, Jan-Hendrik Prinz, Frank Noé, John D. Chodera, and
   Peter G. Bolhuis. "OpenPathSampling: A flexible, open framework for path
   sampling simulations. 2. Building and Customizing Path Ensembles and
   Sample Schemes." J. Chem. Theory Comput. **15**, 837 (2019).
   https://doi.org/10.1021/acs.jctc.8b00627

Both citations: :download:`openpathsampling.bib
<openpathsampling.bib>` (citation keys ``ops1`` and ``ops2``).

--------------------------------------------------------------------------------

For Users
---------

.. toctree::
    :maxdepth: 2

    install
    examples/index
    topics/index
    videos
    cli/index
    ecosystem
    faq
    getting_help


For Developers
--------------

.. toctree::
    :maxdepth: 2
    
    developers/index
    api_sections

.. toctree::
    :hidden:

    acknowledgments

--------------------------------------------------------------------------------


License
-------
OpenPathSampling is licensed under the `MIT license
<https://github.com/openpathsampling/openpathsampling/blob/master/LICENSE>`_.
