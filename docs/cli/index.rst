.. _cli:

Command Line Interface
======================

A separate command line tool for OpenPathSamplng can be installed. It is
available via either ``conda`` or ``pip``:

.. code:: bash

    conda install -c conda-forge openpathsampling-cli
    # or
    python -m pip install openpathsampling-cli

Once you install this, you'll have access to the command
``openpathsampling`` in your shell (although we recommend aliasing that to
either ``paths`` or ``ops`` -- save yourself some typing!)

This command is a gateway to many subcommands, just like ``conda`` and
``pip`` (which have subcommands such as ``install``) or ``git`` (which has
subcommands such as ``clone`` or ``commit``). You can get a full listing all
the subcommands with ``openpathsampling --help``. For more information on
any given subcommand, use ``openpathsampling SUBCOMMAND --help``, replacing
``SUBCOMMAND`` with the subcommand you're interested in.


.. toctree::
    :titlesonly:
    :maxdepth: 1

    workflow
    command_usage
    commands
    compile/index
    expression_eval
    
