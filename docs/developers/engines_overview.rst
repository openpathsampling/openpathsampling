Overview of Engines
===================

From the OpenPathSampling perspective, an engine must involve at least two
parts: it must have a subclass of :class:`.Snapshot`, which describes the
state of the system at a give point in time, and it must have a subclass of
:class:`.DynamicsEngine`, which propagates the system through time. It also
must create a :class:`.SnapshotDescriptor`, which tells the storage
subsystem how large each snapshot will be (e.g., how many atoms).

In practice, most engines have several other objects which help the
:class:`.DynamicsEngine`. Frequently there are some objects which help
describe the interactions: in molecular mechanics, these are "force fields"
and "topologies," or sometimes just the "potential energy surface." There is
usually some customizable rule for propagating the system through time,
typically called an "integrator."

OPS is indifferent to the details of how these sorts of things are
implemented in any specific engine. What matters is that the
:class:`.DynamicsEngine` obeys the OPS API, and that OPS can interact with
the :class:`Snapshots <.Snapshot>`.

The two engine APIs
-------------------

OpenPathSampling supports two approaches to controlling an engine: through
direct control (where a Python API enables you to directly control the MD
and run until the next saved snapshot) and through indirect control (where
the engine runs as an external process).

In general, if your engine can be implemented using the direct control API,
we would suggest using it. However, many existing engines are not designed
to be controlled interactively by a Python script. For these, we have the
indirect control API, which uses the output of a trajectory file as an
intermediary between the engine and OPS.

Customizing snapshots
---------------------

The :class:`.Snapshot` includes all the relevant information about the state
of the system at a given instant in time. For molecular dynamics, this
typically includes positions and velocities, as well as box vectors for
periodic systems. In OPS, it is easy to add other information to snapshots,
known as "snapshot features." This can include storing new information in
each snapshot (such as wavefunction data), computing data based on other
information in the snapshot (such as kinetic energy), or accessing other
information (by reference) that doesn't change between snapshots (such as
particle passes). Implementing particular snapshot features will also
immediately enable more advanced techniques in OPS (such as two-way
shooting, or support for interactions with MDTraj). This is designed to help
you, as an engine developer, get the most functionality with the least
effort.
