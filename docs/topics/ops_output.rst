.. ops-output::

=========================
Output during simulations
=========================

*Or: How do I know that my simulation is running?*

OpenPathSampling is designed to be a library, meaning that we intend many
parts of the code to be reused within other packages that our users might
develop. Indeed, *every* OPS simulation is technically an application that
uses OPS as a library. This means that, for the most part, we try to keep
default output minimal. The primary output from an OPS simulation is the
netCDF file, where all data about the simulation is stored. Other output
along the way is secondary.

While your simulation is running, there are four ways you can see evidence
that it is running. This document will discuss each of them, including some
caveats about how you should or should not use them and how to customize
their behavior.


Progress output (default stdout)
--------------------------------

The most obvious evidence that your simulation is running is our default
output to ``stdout``, which, in the case of a :class:`.PathSampling`
simulation,  includes the Monte Carlo step number, the elapsed time, and an
estimate of the time remaining. An example output:

.. code-block:: text

    Working on Monte Carlo cycle number 26
    Running for 4 minutes 15 seconds - 10.21 seconds per step
    Estimated time remaining: 1 day 4.28 hours

If you don't want to print this to ``stdout``, you can change the value of
:attr:`.PathSampling.output_stream` to any file handler you desire. For
example, if you want to silence it, you can use ``sampler.output_stream =
open(os.devnull)`` (after creating ``sampler`` as an instance of
:class:`.PathSampling`). Note that this also true for any subclass of
:class:`.PathSimulator`, not just :class:`.PathSampling`.


Read the netCDF file
--------------------

One of the reasons we use netCDF as the file format for OPS is that it can
be read while it is open for writing. Since the analysis of OPS data uses
this file, this is the recommended approach. You can begin to develop your
specialized analysis tools while you simulation is still running! The netCDF
file is the primary output of an OPS simulation.

Note that there are a few caveats here: while you can simultaneously have an
open file handler for writing and one for reading, you cannot open a file
for reading at the exact moment that it is being written to. This isn't a
problem for large systems (where you may do one quick block of writing per
hour or longer), but can be a problem for small systems.

See various documentation and examples of analysis to see what you can do
with the resulting netCDF file.


Logging what the code is doing
------------------------------

Another way to see what is happing during the simulation is to enable
logging. Logging will give details of what step we're on, what path mover
has been selected, how far the trajectory has progressed (if generating a
trajectory), and whether the move is accepted.

.. note::
    Do *not* use logging output for analysis. Logging output is not
    considered to be part of the API, and the format and details provided
    are subject to change at any time.

The logging facilities use the standard Python :mod:`logging` library. This can
be configured either in code or using a configuration file. One example
configuration file is our `default info-level configuration
<https://github.com/openpathsampling/openpathsampling/blob/master/examples/resources/logging.conf>`_.

To use it, add the following lines in your code when you would like to
enable logging (after putting the ``logging.conf`` file in your working
directory):

.. code:: python

    import logging.config
    logging.config.fileConfig("logging.conf", disable_existing_loggers=False)

By default, that logs info-level information to a file called
``ops_output.log``. You can change the file (including outputting to
``stdout``), or the format of the log entries, or the verbosity of the logs,
all using standard Python tools.

For more details on how to use Python's logging facilities, see:

* https://docs.python.org/3/library/logging.html
* https://docs.python.org/3/howto/logging.html
* https://docs.python.org/3/howto/logging-cookbook.html

In particular, details on how to use a logging configuration file, see:

* https://docs.python.org/3/howto/logging.html#configuring-logging

.. note::

    You can *not* combine sending progress output to ``stdout`` and sending
    logging information to ``stdout``. When written to ``stdout``, the
    progress information tells your terminal to delete and overwrite the
    lines from the preceding update; combined with logging it will delete
    lines from the logging instead!

See a live visualization of the simulation
------------------------------------------

The last way to see what is happening during you simulation is perhaps the
most fun, but also the least practical. You can visualize the last step of
a path sampling simulation by creating a :class:`.StepVisualizer2D` object,
which projects your paths into the plane of an arbitrary pair of collective
variables. The direction of the path is indicated with a dot as the final
frame (like an arrowhead). The color of the path indicates its ensemble.
Heavy-width trajectories show the current state; light width trajectories
(with hollow final frames) indicate rejected trial moves.

To use a :class:`.StepVisualizer2D` during a path sampling simulation,
assign it to the :attr:`.PathSampling.live_visualization` attribute. It will be
updated after every :attr:`.PathSampling.save_frequency` MC steps -- this is
also how frequently the data is sync'd to disk, and how often we run sanity
checks (ensuring that all paths are in the expected ensembles). By default,
this is after every step, but for performance reasons it is much less frequent
for toy models and small systems.

Of course, this visualization is not practical for long-running simulations,
since it requires an interactive environment. However, the same tool can be
used to replay the simulation from a file. The analysis examples demonstrate
how to do this.
