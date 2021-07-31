import openpathsampling as paths
from openpathsampling.engines.openmm import Snapshot
import copy

try:
    from openmmtools import mcmc
    from openmmtools import states
    import openmmtools
except ImportError:
    HAS_OPENMMTOOLS=False
else:
    HAS_OPENMMTOOLS=True

class EngineNotInitializedError(Exception):
    pass

INITIALIZATION_ERROR = EngineNotInitializedError("Set current_snapshot "
                                                 "before using engine.")

def snapshot_from_sampler_state(sampler_state, engine=None):
    """Generate an OPS SnapShot from a OpenMMTools SamplerState

    Parameters
    ----------
    sampler_state: openmmtools.states.SamplerState
        The ``SamplerState`` containing the positions, velocities, and box
        vectors to define the snapshot.
    engine: :class:`.DynamicsEngine`, optional
        The engine to associate with this snapshot; default is ``None``,
        which will make some functionality (e.g, MDTraj integration)
        impossible.
    """
    return Snapshot.construct(
        coordinates=copy.copy(sampler_state.positions),
        velocities=copy.copy(sampler_state.velocities),
        box_vectors=copy.copy(sampler_state.box_vectors),
        engine=engine
    )

class OpenMMToolsMCEngine(paths.engines.DynamicsEngine):
    """Engine for the OpenMMTools MCMC package.

    This takes user-provided ``thermodynamics_state`` and ``move``, and will
    internally create a ``openmmtools.mcmc.MCMCSampler``, which will be used
    to generate "trajectories."

    Parameters
    ----------
    thermodynamic_state : openmmtools.states.thermodynamic_state
        the thermodynamic state to sample
    move : openmmtools.mcmc.MCMCMove
        the Monte Carlo move to use to propagate the system
    options : Dict
        Dictionary with additional options. Allowed options are:

            'n_steps_per_frame': int, default: 10
                the number of individual steps per snapshot
            'n_frames_max': int, default: 5000
                maximum number of frames before a trajectory aborts

    topology : :class:`.Topology`, optional
        Topology information for this sysetm. If an :class:`.MDTrajTopology`
        is provided, this will enable MDTraj integration.
    """
    def __init__(self, thermodynamic_state, move, options, topology=None):
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
        self.topology = topology
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
        self._sampler.run(self.options['n_steps_per_frame'])
        n_acc_after = self._get_n_accepted()
        if n_acc_before != n_acc_after:
            self.current_snapshot = snapshot_from_sampler_state(
                self._sampler.sampler_state,
                engine=self
            )

        return self.current_snapshot

    def to_dict(self):
        serialize = openmmtools.utils.serialize
        return {
            'thermodynamic_state': serialize(self.thermodynamic_state),
            'move': serialize(self.move),
            'options': self.options,
            'topology': self.topology,
        }

    @classmethod
    def from_dict(cls, dct):
        deserialize = openmmtools.utils.deserialize
        dct['thermodynamic_state'] = deserialize(dct['thermodynamic_state'])
        dct['move'] = deserialize(dct['move'])
        return cls(**dct)

    @property
    def mdtraj_topology(self):
        if self.topology is None:
            raise RuntimeError("Can not convert to MDTraj: no topology "
                               "associated with this engine.")
        return self.topology.mdtraj

