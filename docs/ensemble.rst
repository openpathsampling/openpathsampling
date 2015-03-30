.. _ensemble:

.. currentmodule:: openpathsampling.ensemble

Ensemble Functions
==================

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

Volume Ensembles
----------------------
.. autosummary::
    :toctree: api/generated/

	VolumeEnsemble
	AllInEnsemble
	AllOutEnsemble
        PartInEnsemble
	PartAllOutEnsemble
	
Set-based Ensemble combinations
-----------------------------
.. autosummary::
    :toctree: api/generated/

        IntersectionEnsemble
        UnionEnsemble
        SymmetricDifferenceEnsemble
        RelativeComplementEnsemble
        EnsembleCombination
	
length Ensembles
-----------------------------
.. autosummary::
    :toctree: api/generated/

        LengthEnsemble
    
