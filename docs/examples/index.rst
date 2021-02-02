.. _examples:

Examples
========

This page provides a series of examples, tutorials and recipes to help you
get started with OpenPathSampling.

Each subsection of the pages linked below is a notebook.  To open these
notebooks in a "live" session and execute the documentation interactively,
open them using ``jupyter notebook NOTEBOOK_NAME.ipynb``. 

If you installed OPS from source, you can find these in the ``examples``
directory within the root ``openpathsampling`` directory. You can also find
them in our `GitHub repository
<https://github.com/openpathsampling/openpathsampling/tree/master/examples>`_.

.. note:: It's a *lot* more fun to run these examples live than to just read
          them here!

Tutorial
--------

One of the best ways to get familiar with OPS is our tutorial, which has
been used in several workshops and classrooms over the years. The main part
of the tutorial (notebooks numbered 1-3) is essentially the same as the
alanine dipeptide TPS example, and usually takes no more than 90 minutes to
complete. It can be downloaded at the link below or directly launched in
Binder. Read the README file for details on using the tutorial.

* https://github.com/openpathsampling/ops_tutorial |tutorial-binder|

.. |tutorial-binder| image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/openpathsampling/ops_tutorial/HEAD?urlpath=lab

There's also a `YouTube video of the tutorial 
<https://www.youtube.com/watch?v=a46tOWnqjNY>`_ as it was presented at a
workshop in 2017.


Introductory Examples
---------------------

The next examples give the entire process of a path sampling simulation: going
from an initial frame to a set of initial trajectories, performing the path
sampling, and the analyzing the results.

.. toctree::
   :maxdepth: 2

   AD_tps
   mstis

We recommend beginning with those two examples: they cover most of the
essential points of using OPS. We also have several other examples which
show how to use these approaches for specific cases. Here are several of
those:

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

More examples can be found in the `OPS additional examples repository
<https://github.com/openpathsampling/ops_additional_examples>`_. Live
versions of those examples are also available through `binder
<https://mybinder.org/v2/gh/openpathsampling/ops_additional_examples/HEAD>`_.

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
