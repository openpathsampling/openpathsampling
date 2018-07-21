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

* ``n_steps_per_frame``: the number of internal time steps per saved frame
  (snapshot)
* ``n_frames_max``: maximum number of snapshot allowed in a trajectory
  before the simulation stops

You will need to implement these methods:

* ``current_snapshot`` (``@property``): getter and setter for the current
  snapshot
* ``generate_next_frame()``: takes the (previously set) ``current_snapshot``
  and generates the next saved frame

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

Direct Engine API Step-by-Step
------------------------------

*(a.k.a., "So you want to make a direct API engine?")*

Direct API engines are generally straightforward. 

*TODO: Clarify these steps with more detail*

1. Start by determining what ``features`` your subclass of
   `:class:.Snapshot` will require. Define the snapshot class.

2. Create a snapshot descriptor for your snapshot class.

3. Create the class, add options.

4. Write required methods: ``current_snapshot`` stuff and
   ``generate_next_frame()``

5. Write ``to_dict`` and ``from_dict``, if needed.

6. If the optional methods are worth doing, write them too
