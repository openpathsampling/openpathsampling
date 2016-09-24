OpenPathSampling
================

A Python library to facilitate path sampling algorithms.

OpenPathSampling (OPS) currently supports many variants of transition path
sampling (TPS) and transition interface sampling (TIS), as well as other
useful calculations for rare events, such as committor analysis and flux
calculations.

OPS is independent of the underlying molecular dynamics engine, and
currently has support for OpenMM, as well as an internal engine suitable for
2D toy models.

To learn more about what OPS can do, look at out examples_. If you want to
jump right in, take a look at how easy it is to install_!

.. _examples: examples/index.html
.. _install: getting_started.html

OPS is an open-source project, distributed under the LGPL. Join us in
the development process on GitHub_.

.. _GitHub: http://github.com/openpathsampling/openpathsampling

--------------------------------------------------------------------------------

For Users
---------

.. toctree::
    :maxdepth: 2

    getting_started
    examples/index
    guides/index
    whatsnew
    faq
    developers/index


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
