.. _collectivevariable:

.. currentmodule:: openpathsampling.collectivevariables

Collective Variables
====================

Here we document the various classes that can wrap functions to make them
into OPS collective variables.

For an introduction on collective variables, read the user guide topic on
:ref:`creating-cvs`.

Basic CVs
---------
.. autosummary::
   :toctree: api/generated/

   CollectiveVariable
   FunctionCV
   CoordinateFunctionCV
   GeneratorCV
   CoordinateGeneratorCV
   InVolumeCV

Integrating with other packages
-------------------------------
.. autosummary::
   :toctree: api/generated/

   MDTrajFunctionCV
   PyEMMAFeaturizerCV
   MSMBFeaturizerCV
   PLUMEDCV
   PLUMEDInterface

