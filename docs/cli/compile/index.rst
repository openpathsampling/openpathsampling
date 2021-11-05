Compile
=======

The ``compile`` command takes input files and creates OPS object databases
from them. To a large extent, this plays the same role as a transition input
file for a simulation program: you set the values for certain keywords in a
human-readable file, and the program can create a simulation based on that.

However there are a few significant differences that make the ``compile``
command more flexible and more powerful than a traditional input file:

* **The nested input structure is a better representation of the
  simulation.** Most traditional input files use a flattened structure --
  that is, they simply have .
* **The input file format allows re-use of named objects within a single
  file.**
* **Compiled object databases can be re-used.**


.. ifconfig:: HAS_OPS_CLI is True

    The following sections give information on the input file parameters for
    each type of object.

    .. toctree::
       :maxdepth: 1

       input/engines
       input/collective_variables
       input/volumes
       input/networks
       input/move_schemes
       input/move_strategies
       input/*

.. ifconfig:: HAS_OPS_CLI is False

    The OpenPathSampling CLI does not seem to be installed in this
    environment, so we were unable to build the detailed usage information
    for the ``compile`` command.
