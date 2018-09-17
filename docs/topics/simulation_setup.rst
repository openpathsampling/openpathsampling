.. _simulation_setup:

########################
Setting up path sampling
########################

*Or: What OpenPathSampling does not do for you*

While OpenPathSampling does a lot to simplify the process of performing a
path sampling simulation, it doesn't do everything for you. In particular,
setting up a path sampling simulation can still be challenging.

Just as configurational Monte Carlo requires an initial configuration and
a description of the ensemble (e.g., setting a temperature in an :math:`NVT`
simulation), path sampling requires an initial trajectory and a description
of the path ensemble. In both cases, you use Monte Carlo to obtain the
equilibrated ensemble. However, both describing the ensemble and creating an
initial trajectory are more difficult in path sampling than the equivalent
problems in configurational sampling.

This document discusses the general idea of setting up path sampling
simulations. In includes a few references to specifics of OPS, but is mainly
focused on the concepts. The examples section of the website includes many
practical uses of these ideas.

-----

Defining the ensemble
=====================

The nature of the ensemble depends on the specific algorithm being used.
Here, we'll discuss ensembles for TPS and TIS. In OPS, the :ref:`transition
network <transitions-and-networks>` is used to efficiently create *all* the
ensembles to be used in a simulation. For TPS, there is only one ensemble to
be sampled. For TIS, there is one ensemble per interface (plus additional
helper ensembles, such as the minus ensemble). In multiple state systems,
TIS can involve hundreds of ensembles.

Monte Carlo methods (configurational or path sampling) draw samples from a
statistical ensemble. For configurational sampling, the thermodynamic
ensemble to be sampled is defined by the underlying Hamiltonian of the
system (which includes the masses of the particles and the interaction
energies) and the natural variables of the ensemble (e.g., for the
canonical, or :math:`NVT`, ensemble, these are the number of particles, the
volume, and the temperature).  Similarly, a path ensemble will be defined by
the underlying dynamics (including the underlying thermodynamic ensemble)
and additional restrictions (such as state definitions).

Selecting the extra restrictions for path sampling is much more difficult
that selecting the thermodynamic state to sample: whereas the temperature
and density usually come from experimental considerations, the volume of
phase space that defines a stable state is not an experimentally observable
quantity. Typically, defining the state requires performing MD simulations
in each stable state to characterize its bounds, and there is not a unique
appropriate definition when studying a particular problem of interest.

In TPS, only the states need to be defined. In TIS, one must also define the
interfaces.

Defining good states
--------------------

In TPS, the main traits defining the path ensemble are the state
definitions. The essential point is that states should be stable: a
trajectory started from any point in the state should be expected to take a
long time before visiting any other state.

The topic of state definitions has been included in several review articles.
For more, see:

* |TrajBasedRareEvents|_
* |PractPathSampling|_

.. |TrajBasedRareEvents| replace::
    Bolhuis and Dellago. "Trajectory-based rare event simulations." in
    Reviews in Computational Chemistry, **27** (2010)

.. _TrajBasedRareEvents: http://onlinelibrary.wiley.com/doi/10.1002/9780470890905.ch3/summary

.. |PractPathSampling| replace::
    Bolhuis and Dellago. "Practical and conceptual path sampling issues."
    Eur. Phys. J. Spec. Top. **224**, 209 (2015)

.. _PractPathSampling: https://doi.org/10.1140/epjst/e2015-02419-6


In replica exchange TIS, the role of the minus move in decorrelating
trajectories should also be considered. Ideally, a single minus move
trajectory should be just long enough to fully decorrelate (lose all memory)
within the state. This puts restrictions on both the volume of the state as
well as that of the innermost interface.


Defining good interfaces
------------------------

Much of the literature on TIS uses the expression that "interfaces foliate
phase space." In the OPS picture, where interfaces are volumes instead of
surfaces, this means that inner interfaces (closer to the state) must be
fully contained by outer interfaces. In other words, for every volume in the
ordered set of the state, interface 0, interface 1, etc.; any point in phase
space that is inside interface :math:`i` is also inside interface :math:`i+1`.
Importantly, this includes the state: the innermost interface must entirely
encapsulate the state. Note that OPS does not check this; it requires user
attention.

Interface placement in TIS is essential to its efficiency. If interfaces are
too far apart, the statistics in the estimate for
:math:`P_A(\lambda_{i+1}|\lambda_i)` (see :ref:`tis-analysis`) are poor. If
interfaces are too close together, then to span the transition, we'll need
more interfaces. Since each interface represents an ensemble that must be
sampled to convergence, this can badly hurt the efficiency. Typically, we
aim for about 20-30% overlap between successive interfaces.

There are approaches to `iteratively optimize the location of interfaces
<https://doi.org/10.1063/1.3601919>`_.  These have been implemented for OPS,
and will be added to the core code in the future.

-----

Getting an initial trajectory
=============================

Even with configurational sampling, initial conditions are not always
completely trivial. A PDB entry might be missing hydrogens, and probably
needs to be solvated (and energy-minimized). MC simulations of a liquid
might be started on a lattice, and then melted during an equilibration
phase.  However, obtaining initial conditions for path sampling of rare
events is much more difficult. Path sampling makes it possible to generate
many rare trajectories once a first rare trajectory is given as input. But
how can you get that first rare trajectory?

The answer is, "any way that works." Some approaches include running long
MD, or using a ratcheting approach (such as, or similar to, adaptive
multilevel splitting or forward flux sampling) to get natural trajectories
that go further and further toward the products.

Another important idea is that the initial trajectory does not need to be a
valid physical trajectory from the dynamics of this path ensemble. For
example, it can come from a run at a different temperature. A non-physical
path can be equilibrated to a physical path by performing an equilibration
path sampling simulation before performing the production simulation.

Note that a very bad initial trajectory could equilibrate to a nonsensical
region of trajectory space, just like a bad configuration could equilibrate
to a nonsensical metastable state in configurational Monte Carlo. Selecting
a "nearly" physical path is important.

One approach that has shown promise is to use metadynamics to obtain a first
transition. After the first transition, metadynamics usually will not have
significantly altered the potential energy surface near the barrier, so this
trajectory is likely to equilibrate to a good initial trajectory.

Generating with OPS
-------------------

There are a few approaches within OPS for generating initial trajectories.
One approach is to use run with a higher temperature. This is used in the
alanine dipeptide examples, such as the notebooks for :ref:`AD-tps`. This
can work in practical cases for transitions that are not too rare. However,
here the dynamics is often altered dramatically, and so special attention is
needed. One option is to do a quick committor analysis on the high
temperature trajectories, close to the expected transition. A configuration
with a non-zero committor can be used as a valid good initial trajectory.

Another approach is to use the :class:`.FullBootstrapping` approach, which
starts from a snapshot and rachets up through the path ensembles for a given
transition. This will sample each ensemble until the trajectory satisfies
the next ensemble, then it switches to the new ensemble and samples that.
The process continues until the all ensembles in the transition have
trajectories. It is essentially a version of adaptive multilevel splitting
that has been discretized along the progress parameter. This approach sounds
promising, but in reality is very dependent on the quality of the order
parameter. It is very efficient some simple models, which is why we use it
for the examples like the :ref:`toy-mstis` notebooks. In complicated
systems, it may fail.


Loading from a file
-------------------

In addition to creating a trajectory from scratch using OPS, you can also
load a trajectory from an existing file. This enables you to get an
initial trajectory using whatever tool you're already familiar with (e.g.,
metadynamics with PLUMED). For example, to load a Gromacs XTC file:

.. code-block:: python

    from openpathsampling.engines.openmm.tools import ops_load_trajectory
    traj = ops_load_trajectory("trajectory_file.xtc", top="conf.gro")

The ``top`` argument must be specified as an explicit keyword (using ``=``).
It is required by `MDTraj <http://mdtraj.org>`_, which is used internally to
load files.  Since the MDTraj trajectory object does not have velocities,
OPS will default to giving zero velocity to all atoms. You can still
equilibrate the trajectory by doing two-way shooting with thermalized
velocities (see the `two-way shooting example
<https://gitlab.e-cam2020.eu/dwhswenson/ops_additional_examples/blob/master/two_way_shooting.ipynb>`_
from the ``ops_additional_examples`` repository.)

If your input trajectory file also has velocities, you can assign the
correct velocities to the OPS trajectory by using the
:meth:`.trajectory_from_mdtraj` method. For example, with a Gromacs TRR
file:

.. code-block:: python

    import mdtraj as md
    from openpathsampling.engines.openmm.tools import trajectory_from_mdtraj
    mdt = md.load("trajectory_file.trr", top="conf.gro")
    # next two lines create numpy array for velocities
    trr = md.formats.TRRTrajectoryFile("trajectory_file.trr")
    vel = trr._read(n_frames=len(mdt), atom_indices=None, get_velocities=True)[5]
    # combine `mdt` and `vel` in an OPS trajectory
    traj = trajectory_from_mdtraj(mdt, velocities=vel)

Note that the tricks to get the velocities from the file will depend on what
kind of file you're using. In this particular example, we use part of
MDTraj's private API. The important thing is that the result should be a
Numpy array with shape ``(n_frames, n_atoms, n_spatial)`` and with values in
units of nm/ps.
