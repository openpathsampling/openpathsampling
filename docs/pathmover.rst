.. _pathmover:

.. currentmodule:: opentis.pathmover

PathMover Functions
===================

:class:`opentis.PathMover`

    >>> import openpathsampling as paths
    >>> mover = paths.PathMover()

atomic movers
-------------

These are the leaves of a path tree. The atomic moved that can be made

.. autosummary::
    :toctree: api/generated/

    BackwardShootMover
    ForwardShootMover
    ReplicaExchangeMover


special moves
-------------
.. autosummary::
    :toctree: api/generated/

    MinusMove
    RandomChoiceMover
    PathMover
    PathReversalMover
    EnsembleHopMover
    ReplicaIDChange


collective variable-based volumes
---------------------------------
.. autosummary::
    :toctree: api/generated/

    ConditionalSequentialMover
    PartialAcceptanceSequentialMover
    SequentialMover


mover factory
-------------
.. autosummary::
    :toctree: api/generated/

    PathMoverFactory.NearestNeighborRepExSet
    PathMoverFactory.OneWayShootingSet
    PathMoverFactory.TwoWayShootingSet
