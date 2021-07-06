import abc
import sys
import logging
from future.utils import with_metaclass

import openpathsampling as paths
from openpathsampling.netcdfplus import StorableNamedObject, StorableObject

from ..ops_logging import initialization_logging


logger = logging.getLogger(__name__)
init_log = logging.getLogger('openpathsampling.initialization')

class MCStep(StorableObject):
    """
    A monte-carlo step in the main PathSimulation loop

    It references all objects created and used in a MC step. The used mover,
    and simulator as well as the initial and final sampleset, the step
    number and the generated movechange.

    Attributes
    ----------
    simulation : PathSimulation
        the running pathsimulation responsible for generating the step
    mccycle : int
        the step number counting from the root sampleset
    previous : SampleSet
        the initial (pre) sampleset
    active : SampleSet
        the final (post) sampleset
    change : MoveChange
        the movechange describing the transition from pre to post
    """
    def __init__(self, simulation=None, mccycle=-1, previous=None,
                 active=None, change=None):
        super(MCStep, self).__init__()
        self.simulation = simulation
        self.previous = previous
        self.active = active
        self.change = change
        self.mccycle = mccycle


class PathSimulator(with_metaclass(abc.ABCMeta, StorableNamedObject)):
    """Abstract class for the "main" function of a simulation.

    Parameters
    ----------
    storage : :class:`.Storage`
        Storage file for results

    Attributes
    ----------
    save_frequency : int
        Results should be sync'd (saved to disk) after every
        ``save_frequency`` steps. Note: subclasses must directly implement
        this, the attribute is just a placeholder.
    output_stream : file
        Subclasses should write output to this, allowing a standard way to
        redirect any output.
    allow_refresh : bool
        Whether to allow the output to refresh an ipynb cell; default True.
        This is likely to be overridden when a pathsimulator is wrapped in
        another simulation.
    """
    calc_name = "PathSimulator"
    _excluded_attr = ['sample_set', 'step', 'save_frequency',
                      'output_stream']
    hook_names = ['before_simulation', 'before_step', 'after_step',
                  'after_simulation']

    def __init__(self, storage):
        super(PathSimulator, self).__init__()
        self.storage = storage
        # self.engine = engine
        self.save_frequency = 1
        self.step = 0
        initialization_logging(
            logger=init_log, obj=self,
            entries=['storage']#, 'engine']
        )
        self.sample_set = None
        self.output_stream = sys.stdout  # user can change to file handler
        self.allow_refresh = True
        self.hooks = self.empty_hooks()
        self.attach_default_hooks()

    def sync_storage(self):
        """
        Will sync all collective variables and the storage to disk
        """
        if self.storage is not None:
            self.storage.sync_all()

    def attach_default_hooks(self):
        pass

    def empty_hooks(self):
        """Return a hook dictionary with no hooks."""
        return {k: [] for k in self.hook_names}

    def attach_hook(self, hook, hook_for=None):
        """Attach a hook class or method to this simulation.

        Parameters
        ----------
        hook : :class:`.PathSimulatorHook` or method
            Hook to add
        hook_for : str or None
            If None (default) then the ``hook`` must be a class with methods
            named to match the hook names in this simulator. If ``hook`` is
            a method, then ``hook_for`` must be the name of the hook it
            represents
        """
        def add_hook_method(hook_method, hook_name):
            try:
                self.hooks[hook_name].append(hook_method)
            except KeyError:
                raise TypeError("No hook '" + hook_name + "' in " +
                                str(self.__class__.__name__))

        if hook_for is None:
            for hook_name in hook.implemented_for:
                hook_method = getattr(hook, hook_name)
                add_hook_method(hook_method, hook_name)
        else:
            add_hook_method(hook, hook_for)

    def run_hooks(self, hook_name, **kwargs):
        """Run the hooks for the given ``hook_name``"""
        hook_name_state = {}
        for hook in self.hooks[hook_name]:
            result = hook(**kwargs)
            if result is not None:
                hook_name_state[hook] = result
        if hook_name_state != {}:
            return hook_name_state

    @abc.abstractmethod
    def run(self, n_steps):
        """
        Run the simulator for a number of steps

        Parameters
        ----------
        n_steps : int
            number of step to be run
        """
        raise NotImplementedError()

    def save_initial_step(self):
        """
        Save the initial state as an MCStep to the storage
        """
        mcstep = MCStep(
            simulation=self,
            mccycle=self.step,
            active=self.sample_set,
            change=paths.AcceptedSampleMoveChange(self.sample_set.samples)
        )

        if self.storage is not None:
            self.storage.steps.save(mcstep)
            self.storage.sync_all()


