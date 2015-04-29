.. _volume:

.. currentmodule:: openpathsampling.volume

Volume Functions
================

:class:`openpathsampling.Volume`

    >>> import openpathsampling as paths
    >>> volume = paths.Vnsemble()

basic volumes
-------------
.. autosummary::
    :toctree: api/generated/

    Volume

set-based volume combinations
-----------------------------
.. autosummary::
    :toctree: api/generated/

    FullVolume
    EmptyVolume
    IntersectionVolume
    UnionVolume
    SymmetricDifferenceVolume
    RelativeComplementVolume
    VolumeCombination

collective variable-based volumes
---------------------------------
.. autosummary::
    :toctree: api/generated/

    LambdaVolume
    LambdaVolumePeriodic

voronoi cell volumes
--------------------
.. autosummary::
    :toctree: api/generated/

    VoronoiVolume


volume factory
--------------
.. autosummary::
    :toctree: api/generated/

    VolumeFactory