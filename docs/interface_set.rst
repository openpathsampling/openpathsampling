.. _interface_set:

.. currentmodule:: openpathsampling.high_level.interface_set

Interface Sets
==============

In transition interface sampling and related methods, we need to define a
set of interfaces for each sampling transition. In OPS, this is essentially
just a list of volumes. However, it can be useful for that list to carry
some addition information (such as the associated collective variable), so
the :class:`.InterfaceSet` is a way to package that information.

For TIS, you will usually use the :class:`.VolumeInterfaceSet` or
:class:`.PeriodicVolumeInterfaceSet`.

Abstract classes
----------------
.. autosummary::
   :toctree: api/generated

   InterfaceSet
   GenericVolumeInterfaceSet

Interface Sets for TIS
----------------------
.. autosummary::
   :toctree: api/generated

   VolumeInterfaceSet
   PeriodicVolumeInterfaceSet

