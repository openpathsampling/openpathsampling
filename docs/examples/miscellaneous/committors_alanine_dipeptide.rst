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

The second notebook runs the committor simulation, and the third notebook
analyzes it.

-----

.. notebook:: examples/misc/alanine_dipeptide_committor/1_select_snapshots.ipynb
   :skip_exceptions:

-----

.. notebook:: examples/misc/alanine_dipeptide_committor/2_committor_simulation.ipynb
   :skip_exceptions:

-----

.. notebook:: examples/misc/alanine_dipeptide_committor/3_committor_analysis.ipynb
   :skip_exceptions:

