OpenPathSampling
================

A Python library to facilitate path sampling algorithms. 

OpenPathSampling (OPS) makes it easy to perform many variants of transition
path sampling (TPS) and transition interface sampling (TIS), as well as
other useful calculations for rare events, such as committor analysis and
flux calculations. In addition, it is a powerful library to build new path
sampling methods.

OPS is independent of the underlying molecular dynamics engine, and
currently has support for OpenMM, as well as an internal engine suitable for
2D toy models.

To learn more about what OPS can do, look at our :ref:`examples <examples>`.
If you want to jump right in, take a look at how easy it is to
:ref:`install <install>`!

OPS is an open-source project, distributed under the LGPL. Join us in
the development process on GitHub_, and follow `@pathsampling
<http://twitter.com/pathsampling>`_ on Twitter for updates!

.. _GitHub: http://github.com/openpathsampling/openpathsampling

.. note:: **Project status:** Currently we're at version 0.9.6.  There will
          be some API cleanup before 1.0, but files created with the current
          version will load in 1.0.  Note, however, that files generated
          with Python 2 will not fully load with Python 3 or vice versa
          (partial Python 3 support starting in 0.9.3).

          Overall, OPS is ready for production work. We're already writing
          several papers that used it, in addition to the paper in prep
          about the code itself. Version 1.0 will be released once the paper
          about the code has been accepted for publication. Documentation is
          still in progress; if you have questions, please contact the
          development team.

To see the most recent updates to the code, see the `release notes
<https://github.com/openpathsampling/openpathsampling/releases>`_ page on
GitHub.

--------------------------------------------------------------------------------

For Users
---------

.. toctree::
    :maxdepth: 2

    install
    examples/index
    topics/index
    faq


For Developers
--------------

.. toctree::
    :maxdepth: 2
    
    developers/index
    api_sections

--------------------------------------------------------------------------------


License
-------
OpenPathSampling is licensed under the LGPL, v. 2.1 or later.
