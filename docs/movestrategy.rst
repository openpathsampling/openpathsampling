.. _movestrategy:

.. currentmodule:: openpathsampling.high_level.move_strategy

MoveStrategy API
================

TODO: explain the ideas here

Abstract class
--------------

.. autosummary::
    :toctree: api/generated/

    MoveStrategy


Shooting strategies
-------------------

.. autosummary::
    :toctree: api/generated/

    OneWayShootingStrategy


Replica exchange strategies
---------------------------

.. autosummary::
    :toctree: api/generated/

    NearestNeighborRepExStrategy
    NthNearestNeighborRepExStrategy
    AllSetRepExStrategy
    SelectedPairsRepExStrategy

Replica motion type strategies
------------------------------

.. autosummary::
    :toctree: api/generated/

    ReplicaExchangeStrategy
    EnsembleHopStrategy

Path reversal strategies
------------------------

.. autosummary::
    :toctree: api/generated/

    PathReversalStrategy

Minus move strategies
---------------------

.. autosummary::
    :toctree: api/generated/

    MinusMoveStrategy
    SingleReplicaMinusMoveStrategy

Global organization strategies
------------------------------

.. autosummary::
    :toctree: api/generated/

    OrganizeByMoveGroupStrategy
    OrganizeByEnsembleStrategy
    PoorSingleReplicaStrategy
