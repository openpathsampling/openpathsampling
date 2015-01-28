.. _ensemble:

.. currentmodule:: opentis.ensemble

Ensemble Functions
==================

Trajectory analysis is the heart of MDTraj. These functions can be used to run
a variety of analyses on :class:`opentis.Ensemble` objects.
It's usually as simple as ::

    >>> import opentis as ops
    >>> t = ops.Ensemble()

volume Ensembles
----------------------
.. autosummary::
    :toctree: api/generated/

	EnsembleVolume
	InXEnsemble
	OutXEnsemble
    HitXEnsemble
	LeaveXEnsemble
	
set-based Ensemble combinations
-----------------------------
.. autosummary::
    :toctree: api/generated/

    AndEnsemble
	OrEnsemble
	XorEnsemble
	SubEnsemble
	EnsembleCombination
	
length Ensembles
-----------------------------
.. autosummary::
    :toctree: api/generated/

    LengthEnsemble
    