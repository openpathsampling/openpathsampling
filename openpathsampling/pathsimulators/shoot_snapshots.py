import logging
import openpathsampling as paths

logger = logging.getLogger(__name__)
from .path_simulator import PathSimulator, MCStep

class ShootFromSnapshotsSimulation(PathSimulator):
    """
    Generic class for shooting from a set of snapshots.

    This mainly serves as a base class for other simulation types
    (committor, reactive flux, etc.) All of these take initial snapshots
    from within some defined volume, modify the velocities in some way, and
    run the dynamics until some ensemble tells them to stop.

    While this is usually subclassed, it isn't technically abstract, so a
    user can create a simulation of this sort on-the-fly for some weird
    ensembles.

    Parameters
    ----------
    storage : :class:`.Storage`
        the file to store simulations in
    engine : :class:`.DynamicsEngine`
        the dynamics engine to use to run the simulation
    starting_volume : :class:`.Volume`
        volume initial frames must be inside of
    forward_ensemble : :class:`.Ensemble`
        ensemble for shots in the forward direction
    backward_ensemble : :class:`.Ensemble`
        ensemble for shots in the backward direction
    randomizer : :class:`.SnapshotModifier`
        the method used to modify the input snapshot before each shot
    initial_snapshots : list of :class:`.Snapshot`
        initial snapshots to use
    """
    def __init__(self, storage, engine, starting_volume, forward_ensemble,
                 backward_ensemble, randomizer, initial_snapshots):
        super(ShootFromSnapshotsSimulation, self).__init__(storage)
        self.engine = engine
        # FIXME: this next line seems weird; but tests fail without it
        paths.EngineMover.default_engine = engine
        try:
            initial_snapshots = list(initial_snapshots)
        except TypeError:
            initial_snapshots = [initial_snapshots]
        self.initial_snapshots = initial_snapshots
        self.randomizer = randomizer

        self.starting_ensemble = (
            paths.AllInXEnsemble(starting_volume) & paths.LengthEnsemble(1)
        )

        self.forward_ensemble = forward_ensemble
        self.backward_ensemble = backward_ensemble

        self.forward_mover = paths.ForwardExtendMover(
            ensemble=self.starting_ensemble,
            target_ensemble=self.forward_ensemble
        )
        self.backward_mover = paths.BackwardExtendMover(
            ensemble=self.starting_ensemble,
            target_ensemble=self.backward_ensemble
        )

        # subclasses will often override this
        self.mover = paths.RandomChoiceMover([self.forward_mover,
                                              self.backward_mover])

    def to_dict(self):
        dct = {
            'engine': self.engine,
            'initial_snapshots': self.initial_snapshots,
            'randomizer': self.randomizer,
            'starting_ensemble': self.starting_ensemble,
            'forward_ensemble': self.forward_ensemble,
            'backward_ensemble': self.backward_ensemble,
            'mover': self.mover
        }
        return dct

    @classmethod
    def from_dict(cls, dct):
        obj = cls.__new__(cls)
        # user must manually set a storage!
        super(ShootFromSnapshotsSimulation, obj).__init__(storage=None)
        obj.engine = dct['engine']
        obj.initial_snapshots = dct['initial_snapshots']
        obj.randomizer = dct['randomizer']
        obj.starting_ensemble = dct['starting_ensemble']
        obj.forward_ensemble = dct['forward_ensemble']
        obj.backward_ensemble = dct['backward_ensemble']
        obj.mover = dct['mover']
        return obj


    def run(self, n_per_snapshot, as_chain=False):
        """Run the simulation.

        Parameters
        ----------
        n_per_snapshot : int
            number of shots per snapshot
        as_chain : bool
            if as_chain is False (default), then the input to the modifier
            is always the original snapshot. If as_chain is True, then the
            input to the modifier is the previous (modified) snapshot.
            Useful for modifications that can't cover the whole range from a
            given snapshot.
        """
        self.step = 0
        snap_num = 0
        self.output_stream.write("\n")
        for snapshot in self.initial_snapshots:
            start_snap = snapshot
            # do what we need to get the snapshot set up
            for step in range(n_per_snapshot):
                paths.tools.refresh_output(
                    "Working on snapshot %d / %d; shot %d / %d\n" % (
                        snap_num+1, len(self.initial_snapshots),
                        step+1, n_per_snapshot
                    ),
                    output_stream=self.output_stream,
                    refresh=self.allow_refresh
                )

                if as_chain:
                    start_snap = self.randomizer(start_snap)
                else:
                    start_snap = self.randomizer(snapshot)

                sample_set = paths.SampleSet([
                    paths.Sample(replica=0,
                                 trajectory=paths.Trajectory([start_snap]),
                                 ensemble=self.starting_ensemble)
                ])
                sample_set.sanity_check()
                new_pmc = self.mover.move(sample_set)
                samples = new_pmc.results
                new_sample_set = sample_set.apply_samples(samples)

                mcstep = MCStep(
                    simulation=self,
                    mccycle=self.step,
                    previous=sample_set,
                    active=new_sample_set,
                    change=new_pmc
                )

                if self.storage is not None:
                    self.storage.steps.save(mcstep)
                    if self.step % self.save_frequency == 0:
                        self.sync_storage()

                self.step += 1
            snap_num += 1



class CommittorSimulation(ShootFromSnapshotsSimulation):
    """Committor simulations. What state do you hit from a given snapshot?

    Parameters
    ----------
    storage : :class:`.Storage`
        the file to store simulations in
    engine : :class:`.DynamicsEngine`
        the dynamics engine to use to run the simulation
    states : list of :class:`.Volume`
        the volumes representing the stable states
    randomizer : :class:`.SnapshotModifier`
        the method used to modify the input snapshot before each shot
    initial_snapshots : list of :class:`.Snapshot`
        initial snapshots to use
    direction : int or None
        if direction > 0, only forward shooting is used, if direction < 0,
        only backward, and if direction is None, mix of forward and
        backward. Useful if using no modification on the randomizer.
    """
    def __init__(self, storage, engine=None, states=None, randomizer=None,
                 initial_snapshots=None, direction=None):
        all_state_volume = paths.join_volumes(states)
        no_state_volume = ~all_state_volume
        # shoot forward until we hit a state
        forward_ensemble = paths.SequentialEnsemble([
            paths.AllOutXEnsemble(all_state_volume),
            paths.AllInXEnsemble(all_state_volume) & paths.LengthEnsemble(1)
        ])
        # or shoot backward until we hit a state
        backward_ensemble = paths.SequentialEnsemble([
            paths.AllInXEnsemble(all_state_volume) & paths.LengthEnsemble(1),
            paths.AllOutXEnsemble(all_state_volume)
        ])
        super(CommittorSimulation, self).__init__(
            storage=storage,
            engine=engine,
            starting_volume=no_state_volume,
            forward_ensemble=forward_ensemble,
            backward_ensemble=backward_ensemble,
            randomizer=randomizer,
            initial_snapshots=initial_snapshots
        )
        self.states = states
        self.direction = direction

        # override the default self.mover given by the superclass
        if self.direction is None:
            self.mover = paths.RandomChoiceMover([self.forward_mover,
                                                  self.backward_mover])
        elif self.direction > 0:
            self.mover = self.forward_mover
        elif self.direction < 0:
            self.mover = self.backward_mover

    def to_dict(self):
        dct = super(CommittorSimulation, self).to_dict()
        dct['states'] = self.states
        dct['direction'] = self.direction
        return dct

    @classmethod
    def from_dict(cls, dct):
        obj = super(CommittorSimulation, cls).from_dict(dct)
        obj.states = dct['states']
        obj.direction = dct['direction']
        return obj


