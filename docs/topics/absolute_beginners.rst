.. _absolute-beginners:

############################
OpenPathSampling terminology
############################

We assume that most people interested in OpenPathSampling are familiar with
the basic tools of molecular simulation: molecular dynamics and Monte Carlo
sampling. [#MC]_ However, you might be new to path sampling methods.
Furthermore, OpenPathSampling uses a rather novel perspective on path
sampling in general. This introduction is aimed at users who are either new
to path sampling or who want to have a better understanding of the
conceptual framework that OpenPathSampling is built on.

Path sampling techniques apply, obviously, to "paths." We will use the word
"trajectory" interchangeably with "path." Each path is made up of several
"snapshots" or "frames." Again, we often use these words interchangeably.

Collective variables are functions of a snapshot
================================================

For a system with :math:`N` atoms, each snapshot on a trajectory (in 3D
space) consists of :math:`3N` coordinates and :math:`3N` velocities --- not
to mention extra information about the simulation box. We can't think in
terms of what all those variables are doing, so we usually try to map the
full dynamics to a smaller set of the most important variables. These are
called "collective variables," and we use them to describe the behavior of
molecules in ways that we can better understand.

Of course, there are an infinite number of possible collective variable
definitions, and finding the best ones to describe the process of a given
transition is very difficult. When specifically referring to a collective
variable which is used to describe the progress of a transition, we will
sometimes call that an "order parameter" or a "reaction coordinate."
Selecting a good order parameter is a necessary first step for some of the
methods implemented in OpenPathSampling.

Volumes apply to snapshots
==========================

In OpenPathSampling, we have a class of objects called volumes. These
represent some sort of region in phase space, and they tell you whether a
given snapshot is in that region.

The most common kind of volume is what we call a :class:`.CVDefinedVolume`.
These are volumes based on some collective variable, as described above. In
this case, we set a minimum and maximum value of the collective variable.
Since the collective variable takes a snapshot and maps it to a single
number, we can determine whether any snapshot is within that volume.

Ensembles apply to trajectories
===============================

We also have a class of objects called :class:`.Ensemble`. More correctly,
we should call these "path ensembles" (or even more correctly, "path
ensemble indicators.") The most important function this provides is the
indicator function, which tells you whether a given trajectory is in the
path ensemble or not.  In the same way that a volume takes a snapshot from a
trajectory and tells you whether it includes it or not, an ensemble takes a
*whole* trajectory and tells you whether it is included.

Just as configurational Monte Carlo is done by making changes to each
snapshot, Monte Carlo can also be performed by making changes to
trajectories. This is the idea of transition path sampling (TPS) and later
methods, such as transition interface sampling (TIS), which are implemented
in OpenPathSampling.

In configurational Monte Carlo, you need some function of the configuration
(usually related to the energy) to decide what steps are allowed and with
what probability. In path sampling Monte Carlo, you need some function of
the entire trajectory. The simplest function is the ensemble indicator
function. In the same way that a :class:`.Volume` object can tell you
whether a snapshot is in that :class:`.Volume`, an :class:`.Ensemble` object
can tell you whether or not a trajectory is in that :class:`.Ensemble`.

A new development in OpenPathSampling is the concept of the "can-extend"
functions, ``can_append`` and ``can_prepend``. A given ensemble's
``can_append`` function takes a trajectory, and tells you whether there is
any way that trajectory could be a subtrajectory of a trajectory in the
ensemble.  While these existed implicitly in the stopping conditions from
previous implementations, the new idea in OpenPathSampling is that they, as
well as the ensembles, can be combined in predictable ways to form new and
well-behaved ensembles. See the information on ensembles for more
details.

PathMovers move in path space
=============================

Let's further consider the analogy with configurational Monte Carlo. When
doing importance sampling, you need a way to generate a new trial
configuration from the old configuration. In simple Monte Carlo, this is
often just a random move of one atom. There are also more advanced
approaches such as cluster moves.

For path sampling, we need a similar way to move through path space. This is
what path movers do: they create new trial trajectories for our ensembles.
The most important path mover is probably the shooting move, which exists
(and in implemented in OPS) in many variants. 

.. [#MC] 
   If not, we recommend that you start by reading the first four chapters of
   Frenkel and Smit's "Understanding Molecular Simulation", or find another
   text that provides a similar introduction.

