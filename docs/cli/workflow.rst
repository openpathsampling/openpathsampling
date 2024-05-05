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

The primary CLI tool for setting up simulations is currently ``compile``
command.

The ``compile`` command takes an input file in YAML or JSON format, and
creates an OPS object database based on that input. In many ways, this is
similar to traditional input files for simulation programs -- all the
parameters that define your simulation are provided in a single
human-readable text file. However, there is a slight difference in that the
``compile`` command creates an object database, which can be re-used for
multiple simulations. This makes it easy to ensure that the parameters used
in one simulation are identical to the parameters used in another
simulation.

Both the ``compile`` command only creates new-style
"SimStore" (``.db``) database files, not the older NetCDFPlus (``.nc``)
files. NetCDFPlus support will be dropped in OPS 2.0.

For more on the ``compile`` command, see its detailed documentation:

* :ref:`The Compile Command <compile_command>`

Running Simulations
-------------------

Many of the commands used for running simulations are related to setting
preparing for path sampling. For example, the ``visit-all`` command can be
used with a higher-temperature engine to create trajectories that visit all
stable states. Then the ``equilibrate`` command could be used to prepare
those trajectories for actual path sampling.

The ``md`` command can also be useful in preparing for path sampling by
using it to gauge the stability of proposed stable states.

Analysis
--------

Currently, the recommended way to analyze a simulation is to perform the
analysis interactively in the Jupyter notebook. However, plugins for the CLI
are in development that would allow the CLI to perform analysis.

