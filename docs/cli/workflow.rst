.. _cli-workflow:

Workflow with the CLI
=====================

As always, the process of running a simulation is (1) set up the simulation;
(2) run the simulation; (3) analyze the simulation. The CLI is mainly
focused on step 2, although it also has tools that generally help with OPS
files.

Any stage can be performed by writing Python scripts that interact with the
core OPS library. However, the CLI provides several tools to simplify common
usage of these stages.

.. TODO: add workflow image here?


Simulation Setup
----------------

The CLI has two primary tools for setting up simulations: the ``wizard``
command and the ``compile`` command. The ``wizard`` command provides a
friendly interactive guide to setting up path sampling simulations. This is
particularly useful for beginners.

The ``compile`` command takes an input file in YAML or JSON format, and
creates an OPS object database based on that input. In many ways, this is
similar to traditional input files for simulation programs -- all the
parameters that define your simulation are provided in a single
human-readable text file. However, there is a slight difference in that the
``compile`` command creates an object database, which can be re-used for
multiple simulations. This makes it easy to ensure that the parameters used
in one simulation are identical to the parameters used in another
simulation.

Both the ``wizard`` and ``compile`` functions only create new-style
"SimStore" (``.db``) database files, not the older NetCDFPlus (``.nc``)
files. NetCDFPlus support will be dropped in OPS 2.0.

For more on the ``wizard`` and ``compile`` commands, see the detailed
documentation on each:

* :ref:`The Wizard <wizard_command>`
* :ref:`The Compile Command <compile_command>`

Running Simulations
-------------------

Analysis
--------

Currently, the recommended way to analyze a simulation is to perform the
analysis interactively in the Jupyter notebook. However, plugins for the CLI
are in development that would allow the CLI to perform analysis.

