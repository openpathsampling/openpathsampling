.. _ensemble:

.. currentmodule:: openpathsampling.ensemble

Ensembles API
=============

The concept of the path ensemble, the (appropriately weighted) set of
trajectories satisfying particular conditions, is at the heart of path
sampling techniques.

Path ensembles are created as :class:`openpathsampling.Ensemble` objects.
Making one is usually as simple as ::

    >>> import openpathsampling as paths
    >>> ens = paths.Ensemble()

where you choose the right kind of ensemble and give it the right
initialization parameters. What's more, ensembles can be combined using the
logical infix operators `&` (and) and `|` (or).

Abstract class
--------------

.. autosummary::
   :toctree: api/generated/

   Ensemble

Basic Ensembles
---------------
.. autosummary::
   :toctree: api/generated/

   EmptyEnsemble
   FullEnsemble

Volume Ensembles
----------------
.. autosummary::
   :toctree: api/generated/

   VolumeEnsemble
   AllInXEnsemble
   AllOutXEnsemble
   PartInXEnsemble
   PartOutXEnsemble
   EntersXEnsemble
   ExitsXEnsemble

Set-based Ensemble combinations
-------------------------------
.. autosummary::
   :toctree: api/generated/

   EnsembleCombination
   IntersectionEnsemble
   UnionEnsemble


Length specific Ensembles
-------------------------
.. autosummary::
   :toctree: api/generated/

   LengthEnsemble
   SingleFrameEnsemble
   OptionalEnsemble

Trajectory Altering Ensembles
-----------------------------

.. autosummary::
   :toctree: api/generated/

   ReversedTrajectoryEnsemble
   SuffixTrajectoryEnsemble
   PrefixTrajectoryEnsemble
   WrappedEnsemble

Sequential Ensembles
--------------------
.. autosummary::
   :toctree: api/generated/

   SequentialEnsemble

TIS-specific Ensembles
----------------------
.. autosummary::
   :toctree: api/generated/

   TISEnsemble
   MinusInterfaceEnsemble

Ensemble Functions
------------------
.. autosummary::
   :toctree: api/generated/

   join_ensembles
