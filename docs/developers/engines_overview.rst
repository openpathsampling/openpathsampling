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
direct control (where each step ???) and through indirect control (where the
engine runs as an external process).

In general, if your engine can be implemented using the direct control API,
we would suggest using it. However, many existing engines are not designed
to be controlled interactively by a Python script. For these, we have the
indirect control API, which uses the output of a trajectory file as an
intermediary between the engine and OPS.

Customizing snapshots
---------------------


