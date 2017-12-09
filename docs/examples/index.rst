.. _examples:

Examples
========

This page provides a series of examples, tutorials and recipes to help you
get started with ``openpathsampling``.

Each subsection of the pages linked below is a notebook.  To open these
notebooks in a "live" session and execute the documentation interactively,
open them using ``jupyter notebook NOTEBOOK_NAME.ipynb``. 

If you installed `openpathsampling` from source, you can find these in the
``examples`` directory within the root ``openpathsampling`` directory. You
cn also find them in our `GitHub repository
<https://github.com/openpathsampling/openpathsampling/tree/master/examples>`_.

.. note:: It's a *lot* more fun to run these examples live than to just read
          them here!

Introductory Examples
---------------------

These examples give the entire process of a path sampling simulation: going
from an initial frame to a set of initial trajectories, performing the path
sampling, and the analyzing the results.

.. toctree::
   :maxdepth: 1

   mstis
   AD_tps

We recommend beginning with those two examples: they cover most of the
essential points of using OPS. Additionally, there's a `tutorial version of
the alanine dipeptide TPS example
<http://gitlab.e-cam2020.eu/dwhswenson/ops_tutorial>`_, if you prefer a
version where *you* have to add the missing pieces to make it run.
We also have several other examples which show how to use these approaches
for specific cases. Here are several of those:

.. toctree::
    :maxdepth: 1

    AD_fixed_len_tps
    AD_mstis
    mistis
    srtis

Note that some of those build off of the earlier examples. If working
through the example notebooks yourself, each notebook in a sequence is
numbered so you know the order to run them.

Advanced Examples
-----------------

The advanced examples demonstrate some of the more specialized uses of OPS. 

.. toctree::
    :maxdepth: 1

    AD_tps_advanced_analysis
    custom_movescheme

Special Topics
--------------

In order to illustrate several features of the code, we have also developed
some "special topics" examples. These usually require that you have already
run one of the introductory examples, but demonstrate some optional behavior
you might find interesting.

.. toctree::
    :maxdepth: 1

    storage
    restarting
    splitting_files
    customizing_visualization
    move_schemes_strategies

In addition, there are a set of examples that are hosted outside the core
OPS repository. These are the OPS additional examples, primarily developed
as part of E-CAM, a European Union Horizon 2020 project. Those examples are
available at https://gitlab.e-cam2020.eu/dwhswenson/ops_additional_examples.

Miscellaneous Examples
----------------------

.. toctree::
    :hidden:

    miscellaneous/index

The examples above provide most of the tools that you might need. However,
to document various other tricks and workflows, we have a page of
:ref:`miscellaneous examples <misc-examples>`.

-----

Let us know if you would like to contribute other example notebooks, or have
any suggestions for how these can be improved.
