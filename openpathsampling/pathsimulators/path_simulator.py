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
    #__metaclass__ = abc.ABCMeta

    calc_name = "PathSimulator"
    _excluded_attr = ['sample_set', 'step', 'save_frequency',
                      'output_stream']

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

    def sync_storage(self):
        """
        Will sync all collective variables and the storage to disk
        """
        if self.storage is not None:
            self.storage.sync_all()

    @abc.abstractmethod
    def run(self, n_steps):
        """
        Run the simulator for a number of steps

        Parameters
        ----------
        n_steps : int
            number of step to be run
        """
        pass

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


