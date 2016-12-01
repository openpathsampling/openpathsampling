.. ad_committor::

Committor Simulation and Analysis for Alanine Dipeptide
=======================================================

This sequence of notebooks provides an example of how to calculate the
committor for a real biomolecular system (alanine dipeptide).

The first notebook deals with selecting the snapshots to use as starting
points for the committor. In this particular case, we didn't care about the
statistics of those points, so we select them directly from ``storage``.
This is because these points were being harvested as part of an approach to
estimate the committor at any given point. If you need the statistics of
those points within the TPS ensemble, e.g., for the transition state
ensemble, then you should select your points differently (by looping over
``storage.steps`` and taking the ``active`` sample set).

The second notebook runs the committor simulation. The third notebook
prepares the analysis data in a way that was needed for this particular
project, where each "shot" of the committor was treated as a separate
"experiment." While this approach is not normally needed, if you already
have data that can easily be put into this format, the fourth notebook shows
how it can be used to construct a :class:`.ShootingPointAnalysis` object.
That notebook also proceeds to demonstrate the kinds of analysis that object
can perform.

-----

.. notebook:: examples/misc/alanine_dipeptide_committor/1_select_snapshots.ipynb
   :skip_exceptions:

-----

.. notebook:: examples/misc/alanine_dipeptide_committor/2_committor_simulation.ipynb
   :skip_exceptions:

-----

.. notebook:: examples/misc/alanine_dipeptide_committor/3_committor_analysis.ipynb
   :skip_exceptions:

-----

.. notebook:: examples/misc/alanine_dipeptide_committor/4_analysis_help.ipynb
   :skip_exceptions:

