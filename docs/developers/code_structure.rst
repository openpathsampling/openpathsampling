Code Structure
==============

OpenPathSampling can be divided into several subsystems that interact with
each other. Usually, new additions to the code will focus on only a small
number of these subsystems -- often just one. Knowing how they all fit
together will help you identify what you should add.

-----

Storing to file (Storage)
-------------------------

Underlying everything in OpenPathSampling is the storage subsystem.  Nearly
all OPS objects inherit from our storable object classes, which makes the
objects themselves storable (enabling things like restart files). While you
shouldn't need to modify the implementation of storage, you will need a
basic understand of what storage can and can't automatically do, and you may
occasionally need to override some methods. All of this is covered in
:ref:`the documentation on storage <dev-storage>`.

Representing simulation data
----------------------------

These are the core representations of simulation data.  These are
:class:`.Snapshot`, :class:`.Trajectory`, :class:`.Sample`,
:class:`.MoveChange`, and :class:`.SampleSet`.

You generally won't need to modify or subclass these (except for
:class:`Snapshots <.Snapshot>`, which some :class:`DynamicsEngines
<.DynamicsEngine>` might need to subclass).  If you're writing new analysis
tools, you'll need to be familiar with how these work. See the section on
analyzing OPS output.  ???TODO: WRITE THAT???

Running Dynamics (Engines)
--------------------------

The ``engines`` subsystem runs the dynamics. If you want to add a new
molecular dynamics engine, this is the subsystem you'll need to become
familiar with.

The main classes you may need to subclass are:

* :class:`.DynamicsEngine`
* :class:`.Snapshot`

In addition, you may need to add custom snapshot ``features``. See :ref:`the
documentation on engines <dev-engines>` in order to learn more about adding
support for new engines.

Interpreting Dynamics (CVs and Volumes)
---------------------------------------

For the dynamics to have any sense, we need to be able to analyze each
frame. The subsystem of the code is mainly focused on frame-by-frame
analysis of a trajectory (i.e., the input is a single :class:`.Snapshot`).

The main classes here are:

* :class:`.CollectiveVariable`, which maps a snapshot to some quantity
  (usually a single float). This could be, e.g., a specific distance of
  interest.
* :class:`.Volume`, which maps a snapshot to a boolean value, indicating
  whether the snapshot is within that phase space volume. This could
  represent a stable state, or the "inside" of an interface in TIS.

The main reason you might subclass these is if you want to create a custom
:class:`.CollectiveVariable` wrapper for some analysis package. We already
have wrappers for `MDTraj <http://mdtraj.org>`_ and `PyEmma
<http://emma-project.org>`_; more wrappers would be welcome! Details on what
to do are including in :ref:`the section on collective variables
<dev-collective-variables>`.

Path Ensembles (Ensembles)
--------------------------

For efficient path sampling, we need to not only run dynamics, but also
to know when to stop the dynamics. In the earlier path sampling methods, all
dynamics were of a fixed length. But more recent approaches have gained
significant efficiency by defining stopping criteria that depend on the
trajectory that has been run so far. These stopping criteria can be derived
from the particle path ensemble being sampled, and for us, this is one of
the primary purposes of the :class:`.Ensemble` object.

In general, we discourage you from trying to subclass :class:`.Ensemble`. In
most cases, it is probably easier to use our existing tools and to create
the ensemble you want from them. ???TODO: add docs about writing custom
ensembles???

If you need to create custom ensembles, it is very likely that you will also
want to create a custom :class:`.TransitionNetwork`. See the discussion
under :ref:`"Higher-Level Tools." <higher-level-general>`

Monte Carlo Moves (PathMovers)
------------------------------

Path sampling is Monte Carlo in path space, so of course we need objects
that perform Monte Carlo moves. These are our :class:`PathMovers
<.PathMover>`. 

Creating custom Monte Carlo moves is one of the common tasks for developers
of new techniques in path sampling, and so we have tried to make it easy.
However, in order to make it easy for other users to mix-and-match your
custom path movers with other path movers, you should also implement a
custom :class:`.MoveStrategy`. Details on both of these are in the
:ref:`documentation on subclassing PathMovers and MoveStrategies
<dev-pathmovers-movestrategies>`.

.. _higher-level-general:

Higher-Level Tools (Networks and MoveStrategies)
------------------------------------------------

In real path sampling simulations, there are often many ensembles and many
path movers. Furthermore, frequently each ensemble needs to know its context
-- how it relates to other ensembles -- in order for its analysis to have
any meaning.

For this reason, we have created some higher-level tools which act as
factories for the ensembles and path movers. If you are interested in doing
path sampling of a custom ensemble, then you probably want to write a custom
:class:`.TransitionNetwork` object. This isn't necessary if your ensemble
is only used temporarily, but if you actually intend to sample and analyze
that ensemble, it should be part of the :class:`.TransitionNetwork`. The
documentation on networks ???move from topics??? explains how to create
custom networks.

As mentioned above, if you're writing a custom path mover, you'll also want
to write a custom :class:`.MoveStrategy`. The :class:`.MoveStrategy` allows
a user to define a desired type of move. For example, using
:class:`.MoveStrategy`, a user could, in one or two lines, choose to have a
different method for selecting the shooting point, or could even create two
sets of shooting moves for the different methods.

The :class:`.MoveStrategy` is then appended to the :class:`.MoveScheme`,
which sorts the strategies into a meaningful order, and then asks each
strategy to create its moves. Details on creating :class:`MoveStrategies
<.MoveStrategy>` are with :ref:`the documentation on custom path movers
<dev-pathmovers-movestrategies>`.

Simulations (PathSimulators)
----------------------------

Most of OPS is designed around path sampling: that is, the Monte Carlo
sampling of path ensembles. However, much of the machinery can be used for
other purposes. If you're interested in using OPS for something other than
path sampling, you'll need to create a new :class:`.PathSimulator` subclass.
The :class:`.PathSimulator` is essentially the "main" function of OPS, and
details on subclassing it can be found in :ref:`the section on
PathSimulators <dev-pathsimulators>`.
