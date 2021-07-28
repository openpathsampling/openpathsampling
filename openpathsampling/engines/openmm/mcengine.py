import openpathsampling as paths
from openpathsampling.engines.openmm import Snapshot

try:
    from openmmtools import mcmc
    from openmmtools import states
except ImportError:
    HAS_OPENMMTOOLS=False
else:
    HAS_OPENMMTOOLS=True

class EngineNotInitializedError(Exception):
    pass

INITIALIZATION_ERROR = EngineNotInitializedError("Set current_snapshot "
                                                 "before using engine.")

def snapshot_from_sampler_state(sampler_state, engine=None):
    return Snapshot.construct(
        coordinates=sampler_state.positions.copy(),
        velocities=sampler_state.velocities.copy(),
        box_vectors=sampler_state.box_vectors.copy(),
        engine=engine
    )

class OpenMMToolsMCEngine(paths.engines.DynamicsEngine):
    """Engine for the MCMC utilities in OpenMMTools
    """
    def __init__(self, thermodynamic_state, move, options):
        # we assume that we have standard n_atoms * 3 dimensions
        descriptor = paths.engines.SnapshotDescriptor.construct(
            Snapshot,
            {'n_atoms': thermodynamic_state.system.getNumParticles(),
             'n_spatial': 3}
        )
        super(OpenMMToolsMCEngine, self).__init__(options, descriptor)
        self.thermodynamic_state = thermodynamic_state
        self.move = move
        self.options = options
        # _sampler is the OpenMMTools MCMCSampler
        self._sampler = None
        self._current_snapshot = None

    @property
    def current_snapshot(self):
        if self._current_snapshot is None:
            raise INITIALIZATION_ERROR

        return self._current_snapshot

    @current_snapshot.setter
    def current_snapshot(self, snapshot):
        # early exit if we're setting with the existing snapshot
        if snapshot == self._current_snapshot:
            return

        sampler_state = states.SamplerState(
            positions=snapshot.coordinates,
            velocities=snapshot.velocities,
            box_vectors=snapshot.box_vectors
        )
        if self._sampler is None:
            self._sampler = mcmc.MCMCSampler(
                thermodynamic_state=self.thermodynamic_state,
                sampler_state=sampler_state,
                move=self.move
            )
        else:
            self._sampler.sampler_state = sampler_state
        self._current_snapshot = snapshot


    def _get_n_accepted(self):
        return sum(dct['n_accepted']
                   for dct in self._sampler.move.statistics)


    def generate_next_frame(self):
        if self._sampler is None:
            raise INITIALIZATION_ERROR

        # We only construct a new snapshot if a trial move has been
        # accepted. This way identical snapshots have the same UUID.
        n_acc_before = self._get_n_accepted()
        self._sampler.run(1)
        n_acc_after = self._get_n_accepted()
        if n_acc_before != n_acc_after:
            self.current_snapshot = snapshot_from_sampler_state(
                self._sampler.sampler_state,
                engine=self
            )

        return self.current_snapshot


