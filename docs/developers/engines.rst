.. _dev-engines:

Engines
=======

One of the most powerful features of OpenPathSampling is that it is
independent of the underlying engine. This means that it can be used for a
wide variety of simulations, including many outside its original domain of
molecular dynamics.

Before diving into the details of engines, let's briefly discuss what an
engine does. Path sampling (and related) approaches are based on studying
groups of trajectories. Each trajectory is a time series of "snapshots" of
the system.  In order to generate a trajectory, we need some algorithm to
propagate from one snapshot to the next. This is what the engine does.

OpenPathSampling was developed in the context of molecular dynamics, and
much of the language we use is based on that. If you are unfamiliar with
molecular dynamics, we recommend chapter 4 of Frenkel & Smit. However, OPS
is designed to work with engines far outside the realm of standard molecular
dynamics. Any time series of some sort of point in a multidimensional space
can be studied with OPS. The sections ??? and ??? will touch on how to
implement engines that do not only deal with standard molecular dynamics.

.. toctree::
    :maxdepth: 1

    engines_overview
    snapshot_features
    engines_direct_api
    engines_indirect_api
    advanced_engines
