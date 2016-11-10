.. _storage:

.. currentmodule:: openpathsampling.storage

Storage Functions
=================

:class:`openpathsampling.storage`

    >>> import openpathsampling as paths
    >>> mover = paths.storage.Storage()

storages
--------

These are the leaves of a path tree. The atomic moved that can be made

.. autosummary::
    :toctree: api/generated/

    Storage
    AnalysisStorage

stores
------
.. autosummary::
    :toctree: api/generated/

    SampleStore
    MoveChangeStore
    MCStepStore
    SampleStore
    SampleSetStore
    TrajectoryStore
    BaseSnapshotStore
    FeatureSnapshotStore
