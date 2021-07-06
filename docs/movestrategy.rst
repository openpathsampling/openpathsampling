.. _movestrategy:

.. currentmodule:: openpathsampling.high_level.move_strategy

MoveStrategy
============

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
    ForwardShootingStrategy
    TwoWayShootingStrategy


Replica exchange strategies
---------------------------

.. autosummary::
    :toctree: api/generated/

    NearestNeighborRepExStrategy
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
