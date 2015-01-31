.. _volume:

.. currentmodule:: opentis.volume

Volume Functions
================

:class:`opentis.Volume`

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
    AndVolume
    OrVolume
    XorVolume
    SubVolume
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