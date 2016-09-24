.. _pathsimulators:

Path Simulators
===============

The "main function" in OPS is contained in a :class:`.PathSimulator` object.
Since the primary goal of OpenPathSampling is to perform Monte Carlo
sampling of path space, the :class:`.PathSampling` variant of
:class:`.PathSimulator` has been thoroughly developed for that case.

However, many approaches are very similar to the Monte Carlo path sampling,
but not exactly the same. If you're interested in implementing such a
method, then you should create a new subclass of :class:`.PathSimulator`.

We have already implemented several such subclasses, and if this is your
goal, we recommend you look at those as examples. In particular,
:class:`.CommittorSimulation`, :class:`.DirectSimulation`, and
:class:`.Bootstrapping` will give you an idea of how to implement these.

In general, you need to initialize your custom :class:`.PathSimulator` with
a :class:`.Storage` and any other initialization options you need. The
actual "main function" is called :method:`.run`. Typically, this takes
a number of steps `n_steps` as an argument.

You automatically get an attribute `save_frequency` which should be used in
your `.run` method such that you save evert `save_frequency` simulation
steps.
