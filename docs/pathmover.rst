.. _pathmover:

.. currentmodule:: openpathsampling.pathmover

PathMover Functions
===================

:class:`openpathsampling.PathMover`

    >>> import openpathsampling as paths
    >>> mover = paths.PathMover()

Abstract class
--------------

.. autosummary::
    :toctree: api/generated/

    PathMover

Changing the trajectory
-----------------------

.. autosummary::
    :toctree: api/generated/

    BackwardShootMover
    ForwardShootMover
    PathReversalMover

Changing the ensemble
---------------------

.. autosummary::
    :toctree: api/generated/

    ReplicaExchangeMover
    EnsembleHopMover

Changing the replica ID
-----------------------

.. autosummary::
    :toctree: api/generated/

    ReplicaIDChange

Combining movers
----------------

.. autosummary::
    :toctree: api/generated/

    RandomChoiceMover
    SequentialMover
    ConditionalSequentialMover
    PartialAcceptanceSequentialMover

Pre-made combined movers
------------------------

.. autosummary::
    :toctree: api/generated/

    MinusMove
    OneWayShootingMover


mover factory
-------------
.. autosummary::
    :toctree: api/generated/

    PathMoverFactory.NearestNeighborRepExSet
    PathMoverFactory.OneWayShootingSet
    PathMoverFactory.TwoWayShootingSet
