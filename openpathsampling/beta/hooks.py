"""
Hooks to change class:`.PathSimulator` behavior.

These hooks group several methods together for use as part of a
:class:`.PathSimulator` ``run`` method. They allow for additional
calculations or output at several points in the simulation.
"""
import time
import openpathsampling as paths
from openpathsampling.netcdfplus import StorableNamedObject


class SimulationNotFoundError(RuntimeError):
    """
    Raised when a hook tries to access its parent simulation before knowing it.
    """
    pass


class PathSimulatorHook(StorableNamedObject):
    """Superclass for PathSimulator hooks.

    This implementation is a do-nothing hook. Subclasses should subclass the
    relevant method in order to add hooks PathSimulator objects.
    """
    implemented_for = ['before_simulation', 'before_step', 'after_step',
                       'after_simulation']

    def before_simulation(self, sim, **kwargs):
        pass  # pragma: no-cover

    def before_step(self, sim, step_number, step_info, state):
        pass  # pragma: no-cover

    def after_step(self, sim, step_number, step_info, state, results,
                   hook_state):
        pass  # pragma: no-cover

    def after_simulation(self, sim, hook_state):
        pass  # pragma: no-cover


class StorageHook(PathSimulatorHook):
    """Standard hook for storage.
    """
    implemented_for = ['before_simulation', 'after_step',
                       'after_simulation']
    def __init__(self, storage=None, frequency=None):
        self.storage = storage
        self.frequency = frequency

    def before_simulation(self, sim, **kwargs):
        if self.storage is None:
            self.storage = sim.storage
        if self.frequency is None:
            self.frequency = sim.save_frequency

    def after_step(self, sim, step_number, step_info, state, results,
                   hook_state):
        if self.storage is not None:
            self.storage.save(results)
            if step_number % self.frequency == 0:
                if sim.sample_set is not None:
                    # some PathSimulators never set their sample_set
                    # but PathSimulator.__init__ sets it to None
                    sim.sample_set.sanity_check()
                self.storage.sync_all()

    def after_simulation(self, sim, hook_state):
        if self.storage is not None:
            sim.storage.sync_all()


class ShootFromSnapshotsOutputHook(PathSimulatorHook):
    """Default (serial) output for ShootFromSnapshotsSimulation objects.

    Updates every time a new snapshot is shot from.

    Parameters
    ----------
    output_stream : stream
        where to write the results; default ``None`` uses the simulation's
        ``output_stream``
    allow_refresh : bool
        whether to allow refresh (see :meth:`.refresh_output`); default
        ``None`` uses the simulation's value
    """
    implemented_for = ['before_simulation', 'before_step']
    def __init__(self, output_stream=None, allow_refresh=None):
        self.output_stream = output_stream
        self.allow_refresh = allow_refresh

    def before_simulation(self, sim, **kwargs):
        if self.output_stream is None:
            self.output_stream = sim.output_stream
        if self.allow_refresh is None:
            self.allow_refresh  = sim.allow_refresh

    def before_step(self, sim, step_number, step_info, state):
        snap_num, n_snapshots, step, n_per_snapshot = step_info
        paths.tools.refresh_output(
            "Working on snapshot %d / %d; shot %d / %d" % (
                snap_num+1, n_snapshots, step+1, n_per_snapshot
            ),
            output_stream=self.output_stream,
            refresh=self.allow_refresh
        )


class LiveVisualizerHook(PathSimulatorHook):
    """
    LiveVisualization using the :class:`openpathsampling.StepVisualizer2D`.

    Updates every `simulation.status_update_frequency` MCSteps, where
    simulation is the `PathSimulator` this hook is attached to.

    NOTE: You will have to set PathSimulator.allow_refresh = False
          Otherwise the LiveVisualization will get refreshed away
          (i.e. deleted) right after creation.

    NOTE: Arguments passed to init take precedence over the corresponding
          parameters of the PathSimulator this hook is attached to. They can
          only be accessed through this hook, e.g. as hook.live_visualizer.

    Parameters
    ----------
    live_visualizer : :class:`openpathsampling.StepVisualizer2D`
        default `None` uses the simulations live_visualizer
    status_update_frequency : int
        number of steps between two refreshs of the visualization;
        default `None` uses the simulations value (PathSampling default=1)
    """
    # NOTE: we visualize after step, because otherwise the 'next' MCstep
    #       would depend on the 'previous' one just for viualization
    #       this deviates from the previous implementation but avoids
    #       having to pass the previous MCstep to before_step hooks
    implemented_for = ['before_simulation', 'after_step']

    def __init__(self, live_visualizer=None, status_update_frequency=None):
        self.live_visualizer = live_visualizer
        self.status_update_frequency = status_update_frequency
        self._simulation = None

    @property
    def live_visualizer(self):
        if self._live_visualizer is not None:
            return self._live_visualizer
        elif self._simulation is not None:
            return self._simulation.live_visualizer
        else:
            # hejung: I think here we should return None, because that is the
            # default and leaving the default value unchanged should not lead
            # to an error :)
            return None

    @live_visualizer.setter
    def live_visualizer(self, val):
        self._live_visualizer = val

    @property
    def status_update_frequency(self):
        if self._status_update_frequency is not None:
            return self._status_update_frequency
        elif self._simulation is not None:
            return self._simulation.status_update_frequency
        else:
            raise SimulationNotFoundError("'status_update_frequency' has not "
                                          + "been set and no hosting "
                                          + "simulation known to get a value."
                                          )

    @status_update_frequency.setter
    def status_update_frequency(self, val):
        self._status_update_frequency = val

    def before_simulation(self, sim):
        self._simulation = sim

    def after_step(self, sim, step_number, step_info, state, results,
                   hook_state):
        if step_number % self.status_update_frequency == 0:
            # do we visualize this step?
            if self.live_visualizer is not None and results is not None:
                # do we visualize at all?
                self.live_visualizer.draw_ipynb(results)


class PathSamplingOutputHook(PathSimulatorHook):
    """
    Default (serial) output for PathSamplingSimulation objects.

    Updates every `PathSampling.status_update_frequency` MCSteps.

    NOTE: Arguments passed to init take precedence over the corresponding
          parameters of the PathSimulator this hook is attached to. They can
          only be accessed through this hook, e.g. as hook.output_stream.

    Parameters
    ----------
    output_stream : stream
        where to write the results; default ``None`` uses the simulation's
        ``output_stream``
    allow_refresh : bool
        whether to allow refresh (see :meth:`.refresh_output`); default
        ``None`` uses the simulation's value
    status_update_frequency : int
        number of steps between two refreshs of the visualization;
        default `None` uses the simulations value (PathSampling default=1)
    """
    implemented_for = ['before_simulation', 'before_step', 'after_simulation']

    def __init__(self, output_stream=None, allow_refresh=None,
                 status_update_frequency=None):
        self.output_stream = output_stream
        self.allow_refresh = allow_refresh
        self.status_update_frequency = status_update_frequency
        self._simulation = None

    @property
    def output_stream(self):
        if self._output_stream is not None:
            return self._output_stream
        elif self._simulation is not None:
            return self._simulation.output_stream
        else:
            raise SimulationNotFoundError("'output_stream' has not "
                                          + "been set and no hosting "
                                          + "simulation known to get a value."
                                          )

    @output_stream.setter
    def output_stream(self, val):
        self._output_stream = val

    @property
    def allow_refresh(self):
        if self._allow_refresh is not None:
            return self._allow_refresh
        elif self._simulation is not None:
            return self._simulation.allow_refresh
        else:
            raise SimulationNotFoundError("'allow_refresh' has not "
                                          + "been set and no hosting "
                                          + "simulation known to get a value."
                                          )

    @allow_refresh.setter
    def allow_refresh(self, val):
        self._allow_refresh = val

    @property
    def status_update_frequency(self):
        if self._status_update_frequency is not None:
            return self._status_update_frequency
        elif self._simulation is not None:
            return self._simulation.status_update_frequency
        else:
            raise SimulationNotFoundError("'status_update_frequency' has not "
                                          + "been set and no hosting "
                                          + "simulation known to get a value."
                                          )

    @status_update_frequency.setter
    def status_update_frequency(self, val):
        self._status_update_frequency = val

    def before_simulation(self, sim):
        self._simulation = sim
        self._initial_time = time.time()

    def before_step(self, sim, step_number, step_info, state):
        nn, n_steps = step_info
        elapsed = time.time() - self._initial_time
        paths.tools.refresh_output(
            "Working on Monte Carlo cycle number " + str(step_number)
            + "\n" + paths.tools.progress_string(nn, n_steps, elapsed),
            refresh=self.allow_refresh,
            output_stream=self.output_stream
        )

    def after_simulation(self, sim):
        paths.tools.refresh_output(
            "DONE! Completed " + str(sim.step) + " Monte Carlo cycles.\n",
            refresh=False,
            output_stream=self.output_stream
        )
