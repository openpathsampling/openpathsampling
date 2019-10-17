Direct API for Engines
======================

The direct API works when your underlying engine has an API to directly
control the dynamics. For example, there's a Python function you can use to
tell the engine to run for 10 time steps, and give you the result. The OPS
direct engine API is then just a wrapper to translate between OPS and the
underlying engine's API.

Direct Engine API Overview
--------------------------

If your engine will use direct control, then your engine class should
inherit from :class:`.DynamicsEngine`. More documentation for the methods
discussed below can be found in the documentation for that class.

Engines take an ``options`` dictionary to provide details of how they run.
All direct control engines need to be able to use the options:

* ``n_steps_per_frame``: Number of internal time steps per saved frame
  (snapshot).
* ``n_frames_max``: Maximum number of snapshots allowed in a trajectory
  before the simulation stops.

You will need to implement these methods:

* ``current_snapshot`` (``@property``): Getter and setter for the current
  snapshot. We assume that the current state of the system (e.g., positions
  and velocities) exists inside your engine; this is the code for
  translating that internal representation to/from the OPS representation as
  a :class:`.Snapshot`.
* ``generate_next_frame()``: Takes the (previously set) ``current_snapshot``
  and generates the next saved frame, i.e., it performs
  ``n_steps_per_frame`` individual time steps.

Optionally, you can override several other methods:

* ``start(snapshot=None)``: Run before each trajectory.
* ``stop(trajectory)``: Run after each trajectory ends.
* ``generate_n_frames(n_frames=1)``: Generate a fixed number of frames at
  once; this can be used for performance in some cases. By default, this
  just calls ``generate_next_frame``, but calls it ``n_frames`` times.
* ``is_valid_snapshot(snapshot)``: Check whether the snapshot is physically
  valid (can be used to check types, look for NaNs, etc.)

Note that you are likely to also need to override the ``to_dict`` and
``from_dict`` methods; see the storage documentation for more details.


How the Direct Engine API Runs
------------------------------

This subsection gives an overview of what order the methods that make up the
direct engine API are called in the course of a trajectory. A simulation can
consist of multiple trajectories, so

* **Start of simulation:** Code that needs to be run before performing any
  dynamics should be part of the standard ``__init__`` of your engine
  subclass.
* **Start of dynamics:** Before a specific dynamics run, we (1) set the
  snapshot with ``engine.current_snapshot = snapshot``, using the current
  snapshot property's setter; and (2) call ``engine.start()``.
* **During dynamics:** During the dynamics, we call
  ``generate_next_frame`` (or possibly ``generate_n_frames``), which will
  typically call the getter for ``engine.current_snapshot`` internally (to
  return the snapshot). The logic within OPS determines whether to ask for
  more frames from the engine.
* **After dynamics:** When the OPS logic says that no more frames are
  required from the engine, we call ``engine.stop(trajectory)``, where
  ``trajectory`` is the trajectory we have created. This is a good place to
  do any per-trajectory clean-up that you need, but may simply be a no-op
  for many direct API engines.


Direct Engine API Step-by-Step
------------------------------

*Or: "So you want to make a direct API engine?"*

Direct API engines are generally straightforward. 

1. Start by determining what ``features`` your subclass of
   :class:`.Snapshot` will require. Define the snapshot class. See
   :ref:`dev-snapshot-features` for details.

2. Create the class, and set up its intialization. You will need to:

   A. Add any options that you use should have some, with appropriate
      default values, in the ``_default_options`` class variable (a 
      dictionary).
   B. Create a :class:`.SnapshotDescriptor` instance for your engine
      instance (and pass it to the :class:`.DynamicsEngine` initialization
      through ``super``.)
3. Write required methods: ``current_snapshot`` (getter and setter) and
   ``generate_next_frame()``.

4. Write ``to_dict()`` and ``from_dict()``, if needed. See documentation on
   storage for more.

5. If the optional methods are needed for your engine, write them.
