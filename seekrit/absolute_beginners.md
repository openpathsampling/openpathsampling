# Conceptual introduction for absolute beginners

We assume that most people interested in OpenPathSampling are familiar with
the basic tools of molecular simulation: molecular dynamics and Monte Carlo
sampling.[^1] However, you might be new to path sampling methods.
Furthermore, OpenPathSampling uses a rather novel perspective on path
sampling in general. This introduction is aimed at users who are either new
to path sampling or who want to have a better understanding of the
conceptual framework that OpenPathSampling is built on.

[^1]: If not, we recommend that you start by reading the first four chapters
of Frenkel and Smit's "Understanding Molecular Simulation", or find another
text that provides a similar introduction.

We assume that you're familiar with the idea that molecular dynamics can be
used to generate trajectories. And we also assume 

Path sampling techniques apply, obviously, to "paths". We will use the word
"trajectory" interchangeably with "path". Each path is made up of several
"snapshots" or "frames". Again, we often use these words interchangeably.

## Collective variables are functions of a snapshot

For a system with $N$ atoms, each snapshot on a trajectory (in 3D space)
consists of $3N$ coordinates and $3N$ velocities --- not to mention extra
information about the simulation box. We can't think in terms of what all
those variables are doing, so we usually try to map the full dynamics to a
smaller set of the most important variables. These are called "collective
variables,"

## Volumes apply to snapshots

In OpenPathSampling, we have a class of objects called volumes. These
represent some sort of region in phase space, and they tell you whether a
given snapshot is in that region.

The most common kind of volume is what we call a `LambdaVolume`. These are
volumes based on some collective variable

## Ensembles apply to trajectories

We also have a class of objects called ensembles. Ensembles 

Just as configurational Monte Carlo is done by making changes to each
snapshot, Monte Carlo can also be performed by making changes to
trajectories. This is the idea of transition path sampling (TPS) and later
methods, such as transition interface sampling (TIS), which are implemented
in OpenPathSampling.

In configurational Monte Carlo, you need some function of the configuration
(usually related to the energy) to decide what steps are allowed and with
what probability. In path sampling Monte Carlo, you need some function of th
entire trajectory. The simplest functions are what are called the ensemble
indicator functions, ???

