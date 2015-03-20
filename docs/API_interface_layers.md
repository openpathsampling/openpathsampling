# Contributor APIs

Perhaps you have a problem which is complicated enough that you need to add
major new functionality to OpenPathSampling. In general, we suggest that you
first see whether the functionality you want to add can be done using the
tools introduced at the power user level, but if not, you can implement
useful subclasses by following these instructions.


## Dynamics Engine API

Many tools already exist to generate the trajectories that our path sampling
methods are based on, and we see no reason to reimplement what has already
been better implemented by others. In addition, the tools we provide in OPS
are independent of the nature of the underlying trajectories.

Because of these to points, we interact with the simulation layer, which
generates trajectories, using a very simple API.

To add support for a new dynamics engine, you must create a subclass of our
`DynamicsEngine` object. This subclass must implement the following features:

* `@property` called `Subclass.current_snapshot`: getting this returns the
  current state of the simulation in the form of an
  `openpathsampling.Snapshot` object. Setting initializes a simulation in
  the state given by the `openpathsampling.Snapshot`.
* `next_frame()`: a function that returns the next saved frame from the
  simulation
* `stop()`: a function that tells the simulation to stop. For simulations
  where the frame generation is under direct control, this is often a no-op.
  However, if your simulation uses an external process to generate frames,
  it is likely that you'll continue generating frames after reaching the
  stopping point. This function is where you do the clean-up.
* The type signature for `__init__` must include arguments `filename` (for
  the netCDF output file), `mode` (which determines whether the file should
  be written or just read), and `opts` (a dictionary which contains all the
  options for the engine). Calling `super` with these arguments will open
  the storage file and will also set all the options in `opts` as attributes
  of the object.

The most obvious way to combine OPS with another engine is if that engine
has a convenient API such that you can, in one function call, obtain the
next frame. This is what we mean by "direct control." For a very simple
example of direct control, see our `ToyDynamicsEngine`. The `OpenMMEngine`
provides a much more powerful (and more complicated) interface with direct
control.

However, in some situations you can't exercise direct control over the
simulation: it's quite possible that your dynamics engine doesn't provide a
Python API. Such situations are more complicated; examples will be developed
in the next version of OPS.

## Ensemble API

Path sampling approaches are fundamentally based on sampling path ensembles.
Before trying to create your own `Ensemble` subclass, we strongly recommend
you see if your goal can be met by creative usage of the path ensemble
classes we've already created. However, if that isn't possible, you can
always create a new subclass.

Such an object must be a subclass of `Ensemble`, and it must implement the
following methods:

* `__call__(trajectory)`: determine whether the given trajectory is in the
  ensemble.
* `can_append(trajectory)`: determine whether appending to the given
  trajectory would make it impossible for further appending to be in the
  ensemble. Another way of thinking about this is to call this the logical
  "not" of the forward stop condition: stop propagating the trajectory
  forward if `can_append` returns `False`.
* `can_prepend(trajectory)`: as with `can_append`, but in the direction of
  decreasing time.
* `__not__()`: return the set-logical "not" of the ensemble subclass. That
  is to say, all trajectories (except the zero-length trajectory) should
  either be in `ensemble_subclass(trajectory)` or
  `ensemble_subclass.__not__(trajectory)`. (While not strictly necessary in
  version 1.0, implementation of `__not__()` is highly recommended for
  future compatibility.)

## PathMover API

As with ensembles, we have very powerful facilities for creating new
`PathMover`s without resorting to creating new subclasses. However, if you
must create a new subclass, here's what you should understand about
`PathMover`s.

First, the fundamental object that goes into each move as input is a
`SampleSet`, which is (at heart) a list of `Sample`s. `Sample`s are ???

`PathMover`s can, in general, be split into two categories: "control"
and "change". "Change" `PathMover`s actually generate new `Sample`s, whereas
"Control" `PathMover`s affect which other `PathMover`s are called. See the
power user discussion of `PathMover`s for more details.

??? something about the MovePaths, and then the .move function ???

## Calculation API


