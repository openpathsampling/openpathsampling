import logging
import numpy as np
import openpathsampling as paths

logger = logging.getLogger(__name__)
from .shoot_snapshots import ShootFromSnapshotsSimulation

class ReactiveFluxSimulation(ShootFromSnapshotsSimulation):
    """ Reactive Flux simulations (effective positive flux).

    Parameters
    ----------
    storage : :class:`.Storage`
        the file to store simulations in
    engine : :class:`.DynamicsEngine`
        the dynamics engine to use to run the simulation
    states : list of :class:`.Volume`
        the volumes representing the stable states, first state A then state B
    randomizer : :class:`.SnapshotModifier`
        the method used to modify the input snapshot before each shot
    initial_snapshots : list of :class:`.Snapshot`
        initial snapshots to use
    rc : :class:`.CollectiveVariable`
        reaction coordinate
    """
    def __init__(self, storage, engine=None, states=None, randomizer=None,
                 initial_snapshots=None, rc=None):

        # state definition
        self.states = states
        state_A = states[0]
        state_B = states[1]

        # get min/max reaction coordinate of initial snapshots
        self.rc = rc
        rc_array = np.array(self.rc(initial_snapshots))
        rc_min = np.nextafter(rc_array.min(), -np.inf)
        rc_max = np.nextafter(rc_array.max(), np.inf)

        # define reaction coordinate region of initial snapshots
        # = starting_volume
        self.dividing_surface = paths.CVDefinedVolume(self.rc, rc_min, rc_max)

        # define volume between state A and the dividing surface (including A)
        self.volume_towards_A = paths.CVDefinedVolume(self.rc, -np.inf, rc_max)

        # shoot backward until we hit A but never cross the dividing surface
        backward_ensemble = paths.SequentialEnsemble([
            paths.AllInXEnsemble(state_A) & paths.LengthEnsemble(1),
            paths.AllInXEnsemble(self.volume_towards_A - state_A)
        ])

        # shoot forward until we hit state B without hitting A first
        # caution: since the mover will consist of backward and forward
        #          shoot in sequence, the starting ensemble for the forward
        #          shoot is the output of the backward shoot, i.e. a
        #          trajectory that runs from A to the dividing surface and
        #          not just a point there.
        forward_ensemble = paths.SequentialEnsemble([
            paths.AllInXEnsemble(state_A) & paths.LengthEnsemble(1),
            paths.AllOutXEnsemble(state_A | state_B),
            paths.AllInXEnsemble(state_B) & paths.LengthEnsemble(1),
        ])

        super(ReactiveFluxSimulation, self).__init__(
            storage=storage,
            engine=engine,
            starting_volume=self.dividing_surface,
            forward_ensemble=forward_ensemble,
            backward_ensemble=backward_ensemble,
            randomizer=randomizer,
            initial_snapshots=initial_snapshots
        )

        # create backward mover (starting from single point)
        self.backward_mover = paths.BackwardExtendMover(
            ensemble=self.starting_ensemble,
            target_ensemble=self.backward_ensemble
        )

        # create forward mover (starting from the backward ensemble)
        self.forward_mover = paths.ForwardExtendMover(
            ensemble=self.backward_ensemble,
            target_ensemble=self.forward_ensemble
        )

        # create mover combining forward and backward shooting,
        # abort if backward mover fails
        self.mover = paths.NonCanonicalConditionalSequentialMover([
            self.backward_mover,
            self.forward_mover
        ])

    def to_dict(self):
        ret_dict = {
            'states' : self.states,
            'dividing_surface' : self.dividing_surface,
            'volume_towards_A' : self.volume_towards_A,
            'rc' : self.rc
        }
        return ret_dict

    @classmethod
    def from_dict(cls, dct):
        rf = cls.__new__(cls)

        # replace automatically created attributes with stored ones
        rf.states = dct['states']
        rf.dividing_surface = dct['dividing_surface']
        rf.volume_towards_A = dct['volume_towards_A']
        rf.rc = dct['rc']
        return rf
