.. _cli:

Command Line Interface
======================

A separate command line tool for OpenPathSamplng can be installed. It is
available via either ``conda`` or ``pip``:

.. code:: bash

    conda install -c conda-forge openpathsampling-cli
    # or
    pip install openpathsampling-cli

Once you install this, you'll have access to the command
``openpathsampling`` in your shell (although we recommend aliasing that to
either ``paths`` or ``ops`` -- save yourself some typing!)

This command is a gateway to many subcommands, just like ``conda`` and
``pip`` (which have subcommands such as ``install``) or ``git`` (which has
subcommands such as ``clone`` or ``commit``). You can get a full listing all
the subcommands with ``openpathsampling --help``. For more information on
any given subcommand, use ``openpathsampling SUBCOMMAND --help``, replacing
``SUBCOMMAND`` with the subcommand you're interested in.

Here, we will provide a description of a few of the subcommands that the CLI
tool provides. This documentation may not be fully up-to-date with the more
recent releases of the CLI, so use the CLI help tools to get a fuller
understanding of what is included.

For more details on how the CLI interprets its arguments, and to learn how
to develop plugins for the CLI, see its documentation.  The CLI subcommands
are defined through a plugin system, which makes it very easy for developers
to create new subcommands.

* CLI documentation: https://openpathsampling-cli.readthedocs.io/
* CLI code repository: https://github.com/openpathsampling/openpathsampling-cli/

Workflow with the CLI
---------------------

As always, the process of running a simulation is (1) set up the simulation;
(2) run the simulation; (3) analyze the simulation. The CLI is mainly
focused on step 2, although it also has tools that generally help with OPS
files.

To use it, you'll want to first set up your simulation, e.g., in a Jupyter
notebook, and save the simulation objects to a storage file with
``storage.save(obj)``. You should tag you initial conditions with, e.g.,
``storage.tags['initial_conditions'] = init_conds``.

Close that storage file when you're ready to run the simulation, and then
you can run it with, e.g.,

.. code:: bash

    openpathsampling pathsampling setup.nc -o output.nc --nsteps 1000

If the file you created, ``setup.nc`` has exactly one move scheme and has a
sample set tagged as ``initial_conditions``, the CLI will figure that out.
Your output will be in ``output.nc``, and you can analyze it in a Jupyter
notebook, just as before.

Note that this can be especially useful when computing remotely. You can set
up on your local machine, and then you just need to transfer the
``setup.nc`` to the remote machine (assuming you use internally-stored
snapshots).

Finding your way around the CLI
-------------------------------

Like many command line tools, the OPS CLI has the options ``-h`` or
``--help`` to get help. If you run ``openpathsampling --help`` you should
see something like this

.. code:: none

    Usage: openpathsampling [OPTIONS] COMMAND [ARGS]...

      OpenPathSampling is a Python library for path sampling simulations. This
      command line tool facilitates common tasks when working with
      OpenPathSampling. To use it, use one of the subcommands below. For
      example, you can get more information about the pathsampling tool with:

          openpathsampling pathsampling --help

    Options:
      --log PATH  logging configuration file
      -h, --help  Show this message and exit.

    Simulation Commands:
      visit-all     Run MD to generate initial trajectories
      equilibrate   Run equilibration for path sampling
      pathsampling  Run any path sampling simulation, including TIS variants

    Miscellaneous Commands:
      contents     list named objects from an OPS .nc file
      append       add objects from INPUT_FILE  to another file

The ``--log`` option takes a logging configuration file (e.g., `logging.conf
<https://github.com/openpathsampling/openpathsampling/blob/master/openpathsampling/resources/logging.conf>`_,
and sets that logging behavior. If you use it, it must come before the
subcommand name.

You can find out more about each subcommand by putting ``--help`` *after*
the subcommand name, e.g., ``openpathsampling pathsampling --help``, which
returns

.. code:: none

    Usage: openpathsampling pathsampling [OPTIONS] INPUT_FILE

      General path sampling, using setup in INPUT_FILE

    Options:
      -o, --output-file PATH  output ncfile  [required]
      -m, --scheme TEXT       identifier for the move scheme
      -t, --init-conds TEXT   identifier for initial conditions (sample set or
                              trajectory)
      -n, --nsteps INTEGER    number of Monte Carlo trials to run
      -h, --help              Show this message and exit.

Here you see the list of the options for the running a path sampling
simulation. In general, path sampling requires an output
file, a move scheme and initial conditions from some input file, and the
number of steps to run.  Note that only the output file is technically
required: the CLI will default to running 0 steps (essentially, testing the
validity of your setup), and it can try to guess the move scheme and initial
conditions.  In general, the way it guesses follows the following path:

1. If there is only one object of the suitable type in the INPUT_FILE, use
   that.
2. If there are multiple objects of the correct type, but only one has a
   name, use the named object.
3. In special cases it looks for specific names, such as
   ``initial_conditions``, and will use those.

Full details on how various CLI parameters search the storage can be seen in
the `Parameter Interpretation
<https://openpathsampling-cli.readthedocs.io/en/latest/interpretation.html>`_
section of the CLI docs.

Simulation Commands
-------------------

One of the main concepts when working with the CLI is that you can create
all the OPS simulation objects without running the simulation, save them in
an OPS storage file, and then load them again to actually run your
simulation. For simulation commands, the options all deal with loading
simulation objects from storage.

Here are some of the simulation commands implemented in the OPS CLI:

* ``visit-all``: create initial trajectories by running MD until all states
  have been visited (works for MSTIS or any 2-state system); must provide
  states, engine, and initial snapshot on command line
* ``equilibrate``: run equilibration for path sampling (until first
  decorrelated trajectory); must provide move scheme and initial conditions
  on the command line
* ``pathsampling``: run path sampling with a given move scheme (suitable for
  custom TPS schemes as well as TIS/RETIS); must provide move scheme,
  iniital conditions,  and number of MC steps on command line

Miscellaneous Commands
----------------------

Even for users who prefer to develop their OPS projects entirely in Python,
foregoing the CLI tools to run simulations, some of the "miscellaneous"
commands are likely to be quite useful. Here are some that are available in
the CLI:

* ``contents``: list all the named objects in an OPS storage, organized by
  store (type); this is extremely useful to get the name of an object to use
* ``append`` : add an object from once OPS storage into another one; this is
  useful for getting everything into a single file before running a
  simulation

Customizing the CLI
-------------------

The OPS CLI uses a flexible plugin system to enable users to easily add
custom functionality. This way, you can create and distribute custom
plugins, giving more functionality to other users who would benefit from it,
without adding everything to the core package and thus overwhelming new
users.

Installing a plugin is easy: just create the directory
``$HOME/.openpathsampling/cli-plugins/``, and copy the plugin Python script
into there. For details on how to write a CLI plugin, see the `CLI
development docs <https://openpathsampling-cli.readthedocs.io/>`_.
