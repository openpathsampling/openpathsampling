.. _dev-pathsimulators:

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
:class:`.FullBootstrapping` will give you an idea of how to implement these.

In general, you need to initialize your custom :class:`.PathSimulator` with
a :class:`.Storage` and any other initialization options you need.
:meth:`.PathSimulator.__init__()` takes the :class:`Storage` as a parameter.
The actual "main function" is called :meth:`run() <.PathSimulator.run>`.
Typically, this takes a number of steps ``n_steps`` as an argument.

By inheriting from :class:`.PathSimulator` (and calling its ``__init__``
method) you automatically get a couple useful attributes:

* :attr:`save_frequency <.PathSimulator.save_frequency>`: Use this to
  control how often your subclass saves to disk
* :attr:`output_stream <.PathSimulator.output_stream>`: Any output from your
  subclass should write to this.
