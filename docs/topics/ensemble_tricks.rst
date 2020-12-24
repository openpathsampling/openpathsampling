
.. _ensemble-tricks:

.. currentmodule:: openpathsampling

==============================
Advanced Uses of OPS Ensembles
==============================

One of the novelties of OpenPathSampling is the new set-theoretic treatment
of path ensembles. The second paper on OPS contains more details. This
section will illustrate a few uses of that approach.


Visit all states
----------------

We want a path which contains at least one frame in each state. The question
is, what ensemble can we use to create such a trajectory?

The first obvious thought would be:

.. code::

    goal_ensemble = PartInXEnsemble(A) & PartInXEnsemble(B)

(which can, of course, be further generalized to more states). However,
while that *is* the ensemble we want to eventually satisfy, we can't use its
``can_append`` to create it, because its ``can_append`` always returns
``True``: the trajectory will go on forever!

But we can use a trick: since what we want is the first trajectory that
satisfies ``goal_ensemble``, we know that every shorter trajectory will not
satisfy it. This means that the shorter trajectories must satisfy the
*complement* of ``goal_ensemble``, and the trajectory we want will be the
first trajectory that does *not* satisfy the complement!

So the trick we'll use is to build the trajectory by using the fact that the
shorter trajectories are in the complement of ``goal_ensemble``, which is
given by:

.. code::

    complement = AllOutXEnsemble(A) | AllOutXEnsemble(B)

The ``generate`` function will stop when that is no longer true, giving us the
trajectory we want. This can be directly generalized to more states.

Here we're not even using the ``can_append`` function. That happens to be
the same as the ensemble itself for this particular ensemble, but
conceptually, we're actually using the test of whether a trajectory is in
the ensemble at all.

.. code::

    init_traj_ensemble = paths.AllOutXEnsemble(A) | paths.AllOutXEnsemble(B)
    trajectory = engine.generate(engine.current_snapshot, [init_traj_ensemble])

Those two lines are the entirety of what you need to do to generate a
trajectory that visits both states.

Note that this functionality has been implemented in OPS as the
:class:`.VisitAllStatesEnsemble`. When using that, you *should* use the
``can_append`` method. That ensemble adds progress reporting functionality
to show how many frames have been run and what states have been found so
far, so it is more pleasant to use. However, this idea is what underlies its
implementation.
