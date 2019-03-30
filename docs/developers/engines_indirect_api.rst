Indirect API for Engines
========================

The indirect API launches a trajectory as a separate, and monitors the
trajectory file it outputs until we hit the stopping criterion. Then we kill
the trajectory and clean up. The main things you need to implement for the
OPS indirect engine API are ways to handle reading and writing files for the
engine.

Indirect Engine API Overview
----------------------------

If you intend to use indirect control, your engine class should inherit from
:class:`.ExternalEngine`, and you should consider overriding the following
methods:

* ``read_frame_from_file(filename, frame_num)``: reads the frame from the
  external engine's file format
* ``write_frame_to_file(filename, snapshot, mode="a")``: writes the frame to
  the external engine's format (used to initiate trajectories)
* ``engine_command()``: returns a string of the command to be called by the
  operating system
* ``set_filenames(number)``: sets the file names for step ``number``
* ``cleanup()``: does any clean-up tasks needed after killing the subprocess
  (e.g., removing temporary files)

Additionally, you may wish to override the following options:

* ``killsig`` (class variable): the signal sent to terminate the process
  (default is ``signal.SIGTERM``).
* ``default_sleep_ms`` (set in ``options``): time the engines sleeps before
  checking again whether a new frame has been written. Note that an
  :class:`.ExternalEngine` will automatically optimize the sleep time until
  you set the option ``auto_optimize_sleep`` to ``False``. 


How the Indirect Engine API Runs
--------------------------------

This subsection gives an overview of what order the methods that make up the
direct engine API are called in the course of a trajectory. Keep in mind
that the indirect API internally implements the direct API (i.e., OPS
methods actually call the direct API). So frequently the indirect API
mirrors the direct API, but with a different method name and some default
behavior.

* **Start of simulation:** Code that needs to be run before performing any
  dynamics should be part of the standard ``__init__`` of your engine
  subclass.
* **Start of dynamics:** Before launching the external engine, the code will
  create an input file with the initial snapshot. Then it will run the
  ``prepare`` method and start a subprocess based on the ``engine_command``
  method.
* **During dynamics:** The output trajectory file is monitored with by
  occasionally reading the file with the ``read_frame_from_file`` method.
* **After dynamics:** After the stopping criterion has been reached, the
  engine will send the ``killsig`` to the external process, and then follow
  the instructions in ``cleanup`` to finish processing of the trajectory.

Indirect Engine API Step-by-Step
--------------------------------

*Or: "So you want to make an indirect API engine?"*

