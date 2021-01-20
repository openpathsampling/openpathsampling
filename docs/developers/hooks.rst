.. _dev_hooks:

Hooks in Path Simulators
========================

.. note:: 
    Hooks are a beta feature. Details of the API may change, and hooks
    may be removed/significantly modified in future major versions.

In order to provide more flexibility, there are certain points where a
developer can modify the functionality of a path simulator. We refer to
these as "hooks." You can create your own hooks by making a subclass of
:class:`.PathSimulatorHook`.  That class has 4 methods (which default to
no-ops) that you can override to add the desired functionality:

* ``before_simulation``: This performs any setup, such as a initialization
  of the hook itself, before the simulation starts. 
* ``before_step``: This runs before each step. It is primarily for progress
  reporting.
* ``after_step`` This runs after each step. It can perform on-the-fly
  analysis of the results from a simulation step, and can adjust the
  behavior of the simulation based on those results.
* ``after_simulation``: This runs at the end of the simulation.

An instance of a subclass of :class:`.PathSimulatorHook` can be attached to
the simulator with ``simulator.attach_hook(hook)``.

Step info
---------

Each simulation type provides a more detailed breakdown of the step
progress, which is provided to both the ``before_step`` and ``after_step``
methods. This is a tuple, the exact contents of which depend on the specific
simulation type.

For a :class:`.ShootFromSnapshotsSimualtion`, such as a
:class:`.CommittorSimulation`, this is the 4-tuple ``(snap, n_snaps, step,
n_steps)`` where ``snap`` is the snapshot number, ``step`` is the step
number within that snapshot, and ``n_snaps`` and ``n_steps`` are the total
number of each, respectively.

Simulation state
----------------

The ``before_step`` and ``after_step`` stages receive a ``state`` variable,
which represents the current state of the simulation. The specific nature of
this again depends on the specific simulation type. For a
:class:`.ShootFromSnapshotsSimulation`, this will be the current initial
snapshot.

Internal state of hooks
-----------------------

If a hook needs to update some sort of internal state over the course of a
simulation, that should be done as part of the ``after_step`` stage. The
``after_step`` method should also return any state that needs to be
preserved. This should be in the form of a dictionary.
