.. _setting_up_sample_sets:

Setting up sample sets
======================

*Or: How do I get the initial conditions?*

Path sampling methods such as TIS involving simultaneously sampling multiple
path ensembles. This means that we need to not only know the trajectory, but
which ensemble it came from. Because of this, OPS uses data objects
:class:`.Sample` to associate a trajectory with an ensemble and
:class:`.SampleSet` to collect multiple :class:`.Sample` instances.

For simulations such as :class:`.PathSampling`, you must provide a
:class:`.SampleSet` as initial conditions for the simulation. This document
deals with several ways to associate trajectories with ensembles, assuming
you've already generated a valid trajectory. See :ref:`get_init_traj` for
details on generating the trajectory.

The following sections provide several options for how to get a
:class:`.SampleSet` once you've obtained the relevant trajectories and have
the ensemble objects (often contained in a :class:`.TransitionNetwork` that
you've created). We'll discuss the advantages and disadvantages of each
approach.


Loading from a file
-------------------

This is the easiest, and will probably be the one you want to use whenever
possible.

.. code:: python

    storage = paths.Storage("myfile.nc", mode='r')
    # often you'll want to load the state of the last accepted MC step:
    final_step = storage.steps[-1]
    sample_set = final_step.active

    # alternatively, you might want a specific sample set you stored
    sample_set = storage.sample_sets[42]  # if you know you want 42

This returns exactly the sample set that previous existed, including the
connection to the previously-used ensembles.  Although this is probably the
best approach for most use cases, there are important situations where you
would not us it:

1. If you don't already have a file with OPS sample sets (chicken and egg,
   right?)
2. If you don't want to associate the trajectories with the same ensembles as
   before. This might be because you're changing the network that you're
   sampling, e.g., using TPS trajectories as initial conditions for TIS, or
   changing the TIS network you're using.


Using the move scheme
---------------------

The move scheme knows the list of all ensembles that it might require for
the first move, so you should use it to ensure that your sample set includes
representatives for every ensemble. It can also take trajectories and
associate them with ensembles. This is a good approach for creating initial
conditions the first time you set up a simulation.

.. code:: python

    # scheme is a MoveScheme object
    # trajectories is a trajectory or list of trajectories
    sample_set = scheme.initial_conditions_from_trajectories(trajectories)

This will also give some output on missing ensembles/extra ensembles.
Ensembles are considered "missing" if they might be required as input for
the move scheme, but they don't have a trajectory associated with them in
the sample set. Ensembles are considered "extra" if they have a
representative in the sample set, but can't be used by the move scheme (not
possible in this setup process).

Aside: Sanity checks
--------------------

There are a few ways to make sure that your sample set is reasonable for
your simulation. OPS will automatically run these before running a path
sampling simulation, but you can check them yourself. Note that they
function based on ``assert`` statements, so this won't work if you disable
asserts with ``python -O``.

.. code:: python

    # assert that each trajectory can be in the associated ensemble
    sample_set.sanity_check()

    # assert that the sample set has the right ensembles represented to be
    # initial conditions for the move scheme
    scheme.assert_inital_conditions(sample_set)


Other approaches for sample sets
--------------------------------

The first two use cases, loading from a file and using the move scheme's
``initial_conditions_from_trajectories`` method, will probably meet nearly
all of your needs. However, there are a few other approaches. These are
legacy approaches that existed before the more general and simpler
approaches were fully stabilized, but they might still be useful.

Mapping equivalent ensembles
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

All objects in OpenPathSampling have a unique universal identifier (UUID)
that gets set when they are created. However, it is possible to create two
objects (e.g., two ensembles) that are equivalent, but do not share the same
UUID. This would occur if you created the same ensemble in two different
networks (e.g., by creating a new network with fewer ensembles than the
original one).

.. code:: python

    sample_set = paths.SampleSet.translate_ensembles(old_sample_set, new_ensembles)

The main use case where this would make more sense than using the move
scheme would be if you wanted to ensure that the ensembles for each
trajectory was preserved, e.g., continuing a simulation with a modified
network. However, be aware that there's no guarantee that the analysis tools
will correctly handle data that combines results from both networks.

Manually matching trajectories and ensembles
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Of course, you can always manually create samples, and put them into a
sample set:

.. code:: python

    samp0 = paths.Sample(replica=0, trajectory=traj0, ensemble=ens0)
    samp1 = paths.Sample(replica=1, trajectory=traj1, ensemble=ens1)
    ...
    sampN = paths.Sample(replica=N, trajectory=trajN, ensemble=ensN)
    sample_set = paths.SampleSet([samp0, samp1, ..., sampN])

In all cases, we strongly recommend that you double check the correctness of
the sample set using the sanity checks listed above as soon as you've
created the sample sets. This can save later confusion.
