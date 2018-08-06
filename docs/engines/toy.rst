.. _toy:

.. module:: openpathsampling.engines.toy

Toy Engine
==========

The toy engine in OPS is a simple engine written in Python that is primarily
designed to be used on simple 2D toy models.

Main objects
------------

.. autosummary::
   :toctree: ../api/generated/

   ToySnapshot
   Snapshot
   ToyEngine
   Engine
   Topology

Additional snapshot features
----------------------------

.. autosummary::
   :toctree: ../api/generated/

   features.instantaneous_temperature

Integrators
-----------

.. autosummary::
   :toctree: ../api/generated/

   integrators.ToyIntegrator
   LeapfrogVerletIntegrator
   LangevinBAOABIntegrator


Potential energy surface (PES) tools
------------------------------------

.. autosummary::
   :toctree: ../api/generated/

   PES
   PES_Sub
   PES_Add
   HarmonicOscillator
   Gaussian
   OuterWalls
   LinearSlope

