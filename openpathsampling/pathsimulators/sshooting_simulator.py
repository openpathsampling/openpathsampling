import logging
import openpathsampling as paths

logger = logging.getLogger(__name__)
from .shoot_snapshots import ShootFromSnapshotsSimulation

class SShootingSimulation(ShootFromSnapshotsSimulation):
    """ S-Shooting simulations.

    Parameters
    ----------
    storage : :class:`.Storage`
        the file to store simulations in
    engine : :class:`.DynamicsEngine`
        the dynamics engine to use to run the simulation
    state_S : :class:`.Volume`
        the volume representing saddle region S.
    randomizer : :class:`.SnapshotModifier`
        the method used to modify the input snapshot before each shot
    initial_snapshots : list of :class:`.Snapshot`
        initial snapshots to use.
    trajectory_length : int
        Trajectory length l of backward/forward shot, total length of generated
        trajectories is 2*l+1. The harvested subtrajectories have length l+1.
    """
    def __init__(self, storage, engine=None, state_S=None, randomizer=None,
                 initial_snapshots=None, trajectory_length=None):

        # Defintion of state S (A and B are only required for analysis).
        self.state_S = state_S

        # Set forward/backward shot length.
        self.trajectory_length = trajectory_length
        l = self.trajectory_length

        # Define backward ensemble:
        # trajectory starts in S and has fixed length l.
        backward_ensemble = paths.SequentialEnsemble([
            paths.LengthEnsemble(l),
            paths.AllInXEnsemble(state_S) & paths.LengthEnsemble(1)
        ])

        # Define forward ensemble:
        # CAUTION: first trajectory is in backward ensemble,
        # then continues with fixed length l.
        forward_ensemble = paths.SequentialEnsemble([
            paths.LengthEnsemble(l),
            paths.AllInXEnsemble(state_S) & paths.LengthEnsemble(1),
            paths.LengthEnsemble(l)
        ])

        super(SShootingSimulation, self).__init__(
            storage=storage,
            engine=engine,
            starting_volume=state_S,
            forward_ensemble=forward_ensemble,
            backward_ensemble=backward_ensemble,
            randomizer=randomizer,
            initial_snapshots=initial_snapshots
        )

        # Create backward mover (starting from single point).
        self.backward_mover = paths.BackwardExtendMover(
            ensemble=self.starting_ensemble,
            target_ensemble=self.backward_ensemble
        )

        # Create forward mover (starting from the backward ensemble).
        self.forward_mover = paths.ForwardExtendMover(
            ensemble=self.backward_ensemble,
            target_ensemble=self.forward_ensemble
        )

        # Create mover combining forward and backward shooting. No condition
        # here, shots in both directions are executed in any case.
        self.mover = NonCanonicalSequentialMover([
            self.backward_mover,
            self.forward_mover
        ])

    def to_dict(self):
        ret_dict = {
            'state_S' : self.state_S,
            'trajectory_length' : self.trajectory_length
        }
        return ret_dict

    @classmethod
    def from_dict(cls, dct):
        sshooting = cls.__new__(cls)

        # replace automatically created attributes with stored ones
        sshooting.state_S = dct['state_S']
        sshooting.trajectory_length = dct['trajectory_length']
        return sshooting

class NonCanonicalSequentialMover(
          paths.SequentialMover):
    """ Special mover for reactive flux simulation.

    This mover inherits from :class:`.SequentialMover` and
    alters only the `move` method to return the output of the corresponding
    :class:`.NonCanonicalSequentialMoveChange`.
    """
    _is_canonical = False

    def move(self, sample_set):
        change = super(NonCanonicalSequentialMover,
                       self).move(sample_set)
        return NonCanonicalSequentialMoveChange(
            subchanges=change.subchanges,
            mover=change.mover,
            details=change.details
        )

class NonCanonicalSequentialMoveChange(
          paths.SequentialMoveChange):
    """ Special move change for reactive flux simulation.

    This move change inherits from :class:`.SequentialMoveChange`
    and returns the outcome of the last subchange.
    """
    @property
    def canonical(self):
        return self.subchanges[-1]
