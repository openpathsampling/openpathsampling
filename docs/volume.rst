.. _volume:

.. currentmodule:: openpathsampling.volume

Volumes
=======

A :class:`.Volume` in OpenPathSampling is something that defines a
hypervolume in phase space. This can involve coordinates as well as
velocities (or any other snapshot feature, such as box vectors). Since
volumes represent a set of points, they can be combined using set-theoretic
intersections, unions, etc.

Stable states in path sampling are defined in terms of :class:`.Volume`
objects. In addition, OPS defines interfaces in TIS as a :class:`.Volume`.

Basic Volumes
-------------
.. autosummary::
    :toctree: api/generated/

    Volume
    FullVolume
    EmptyVolume

Set-Based Volume Combinations
-----------------------------
.. autosummary::
    :toctree: api/generated/

    IntersectionVolume
    UnionVolume
    SymmetricDifferenceVolume
    RelativeComplementVolume
    VolumeCombination

Collective Variable-Based Volumes
---------------------------------
.. autosummary::
    :toctree: api/generated/

    CVDefinedVolume
    PeriodicCVDefinedVolume
