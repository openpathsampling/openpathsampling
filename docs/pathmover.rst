.. _pathmover:

.. currentmodule:: openpathsampling.pathmover

PathMover API
=============

:class:`openpathsampling.PathMover`

    >>> import openpathsampling as paths
    >>> mover = paths.PathMover()

Abstract class
--------------

.. autosummary::
    :toctree: api/generated/

    PathMover
    SampleMover
    EngineMover
    SelectionMover
    SubtrajectorySelectMover

Changing the trajectory
-----------------------

.. autosummary::
    :toctree: api/generated/

    BackwardShootMover
    ForwardShootMover
    PathReversalMover
    BackwardExtendMover
    OneWayExtendMover
    FinalSubtrajectorySelectMover
    FirstSubtrajectorySelectMover
    RandomSubtrajectorySelectMover

Changing the ensemble
---------------------

.. autosummary::
    :toctree: api/generated/

    ReplicaExchangeMover
    EnsembleHopMover

Combining movers
----------------

.. autosummary::
    :toctree: api/generated/

    RandomChoiceMover
    SequentialMover
    ConditionalSequentialMover
    PartialAcceptanceSequentialMover

Swapping movers
---------------

.. autosummary::
    :toctree: api/generated/

    ReplicaExchangeMover
    StateSwapMover

Logical movers
--------------

.. autosummary::
    :toctree: api/generated/

    ConditionalSequentialMover
    PartialAcceptanceSequentialMover
    SequentialMover
    SubPathMover
    EnsembleFilterMover
    FirstAllowedMover
    LastAllowedMover


Pre-made combined movers
------------------------

.. autosummary::
    :toctree: api/generated/

    MinusMover
    OneWayShootingMover
