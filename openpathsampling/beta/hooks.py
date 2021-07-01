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


def _self_or_sim_property_or_err(obj, prop_name, sim_prop_name=None):
    # helper function to get values either from the hook or the parent sim
    if sim_prop_name is None:
        # make it possible to have different names for the sim and
        # self attributes (but default to same name)
        sim_prop_name = prop_name
    if getattr(obj, "_" + prop_name) is not None:
        return getattr(obj, "_" + prop_name)
    elif obj._simulation is not None:
        return getattr(obj._simulation, sim_prop_name)
    else:
        raise SimulationNotFoundError("'" + prop_name + "' has not "
                                      + "been set and no hosting "
                                      + "simulation known to get "
                                      + "simulation." + sim_prop_name + " ."
                                      )


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
    """
    Standard hook for storage.

    NOTE: Arguments passed to init take precedence over the corresponding
          parameters of the PathSimulator this hook is attached to. They can
          only be accessed through this hook, e.g. as hook.live_visualizer.

    Parameters
    ----------
    storage : :class:`.Storage`
              where to save to; default ``None`` uses the simulation's
              ``storage``
    frequency : int
                save frequency measured in steps; default ``None`` uses the
                simulation's value for ``save_frequency``

    """
    implemented_for = ['before_simulation', 'after_step',
                       'after_simulation']

    def __init__(self, storage=None, frequency=None):
        self.storage = storage
        self.frequency = frequency
        self._simulation = None

    @property
    def frequency(self):
        return _self_or_sim_property_or_err(self, "frequency", "save_frequency")

    @frequency.setter
    def frequency(self, val):
        self._frequency = val

    @property
    def storage(self):
        return _self_or_sim_property_or_err(self, "storage")

    @storage.setter
    def storage(self, val):
        self._storage = val

    def before_simulation(self, sim, **kwargs):
        self._simulation = sim

    def after_step(self, sim, step_number, step_info, state, results,
                   hook_state):
        if self.storage is not None:
            self.storage.save(results)
            if step_number % self.frequency == 0:
                self.storage.sync_all()

    def after_simulation(self, sim, hook_state):
        if self.storage is not None:
            self.storage.sync_all()
        if sim.storage is not None:
            if self.storage is not sim.storage:
                # can/should call sync_all only once per storage
                sim.storage.sync_all()


class ShootFromSnapshotsOutputHook(PathSimulatorHook):
    """Default (serial) output for ShootFromSnapshotsSimulation objects.

    Updates every time a new snapshot is shot from.

    NOTE: Arguments passed to init take precedence over the corresponding
          parameters of the PathSimulator this hook is attached to. They can
          only be accessed through this hook, e.g. as hook.live_visualizer.

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
        self._simulation = None

    @property
    def output_stream(self):
        return _self_or_sim_property_or_err(self, "output_stream")

    @output_stream.setter
    def output_stream(self, val):
        self._output_stream = val

    @property
    def allow_refresh(self):
        return _self_or_sim_property_or_err(self, "allow_refresh")

    @allow_refresh.setter
    def allow_refresh(self, val):
        self._allow_refresh = val

    def before_simulation(self, sim, **kwargs):
        self._simulation = sim

    def before_step(self, sim, step_number, step_info, state):
        snap_num, n_snapshots, step, n_per_snapshot = step_info
        paths.tools.refresh_output(
            "Working on snapshot %d / %d; shot %d / %d" % (
                snap_num+1, n_snapshots, step+1, n_per_snapshot
            ),
            output_stream=self.output_stream,
            refresh=self.allow_refresh
        )


class SampleSetSanityCheckHook(PathSimulatorHook):
    """
    Check sample set sanity.

    Parameters
    ----------
    frequency : int
                check frequency measured in steps; default ``None`` uses the
                simulation's value for ``save_frequency``
    """
    implemented_for = ['before_simulation', 'after_step']

    def __init__(self, frequency=None):
        self.frequency = frequency
        self._simulation = None

    @property
    def frequency(self):
        return _self_or_sim_property_or_err(self, "frequency", "save_frequency")

    @frequency.setter
    def frequency(self, val):
        self._frequency = val

    def before_simulation(self, sim, **kwargs):
        self._simulation = sim

    def after_step(self, sim, step_number, step_info, state, results,
                   hook_state):
        if step_number % self.frequency == 0:
            if sim.sample_set is not None:
                # some PathSimulators never set their sample_set
                # but PathSimulator.__init__ sets it to None
                sim.sample_set.sanity_check()


class LiveVisualizerHook(PathSimulatorHook):
    """
    LiveVisualization using the :class:`openpathsampling.StepVisualizer2D`.

    Updates every ``simulation.status_update_frequency`` MCSteps, where
    simulation is the ``PathSimulator`` this hook is attached to.

    NOTE: You will have to set PathSimulator.allow_refresh = False
          Otherwise the LiveVisualization will get refreshed away
          (i.e. deleted) right after creation.

    NOTE: Arguments passed to init take precedence over the corresponding
          parameters of the PathSimulator this hook is attached to. They can
          only be accessed through this hook, e.g. as hook.live_visualizer.

    Parameters
    ----------
    live_visualizer : :class:`openpathsampling.StepVisualizer2D`
        default ``None`` uses the simulation's live_visualizer
    status_update_frequency : int
        number of steps between two refreshs of the visualization;
        default ``None`` uses the simulation's value (PathSampling default=1)
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
        return _self_or_sim_property_or_err(self, "live_visualizer")

    @live_visualizer.setter
    def live_visualizer(self, val):
        self._live_visualizer = val

    @property
    def status_update_frequency(self):
        return _self_or_sim_property_or_err(self, "status_update_frequency")

    @status_update_frequency.setter
    def status_update_frequency(self, val):
        self._status_update_frequency = val

    def before_simulation(self, sim, **kwargs):
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

    Updates every ``PathSampling.status_update_frequency`` MCSteps.

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
        return _self_or_sim_property_or_err(self, "output_stream")

    @output_stream.setter
    def output_stream(self, val):
        self._output_stream = val

    @property
    def allow_refresh(self):
        return _self_or_sim_property_or_err(self, "allow_refresh")

    @allow_refresh.setter
    def allow_refresh(self, val):
        self._allow_refresh = val

    @property
    def status_update_frequency(self):
        return _self_or_sim_property_or_err(self, "status_update_frequency")

    @status_update_frequency.setter
    def status_update_frequency(self, val):
        self._status_update_frequency = val

    def before_simulation(self, sim, **kwargs):
        self._simulation = sim
        self._initial_time = time.time()

    def before_step(self, sim, step_number, step_info, state):
        if step_number % self.status_update_frequency == 0:
            nn, n_steps = step_info
            elapsed = time.time() - self._initial_time
            paths.tools.refresh_output(
                "Working on Monte Carlo cycle number " + str(step_number)
                + "\n" + paths.tools.progress_string(nn, n_steps, elapsed),
                refresh=self.allow_refresh,
                output_stream=self.output_stream
            )

    def after_simulation(self, sim, hook_state):
        paths.tools.refresh_output(
            "DONE! Completed " + str(sim.step) + " Monte Carlo cycles.\n",
            refresh=False,
            output_stream=self.output_stream
        )
