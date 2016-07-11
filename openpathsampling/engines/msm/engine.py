import numpy as np

from openpathsampling.engines import DynamicsEngine, SnapshotDescriptor
from snapshot import MSMSnapshot as Snapshot


class MSMEngine(DynamicsEngine):
    """
    The trick is that we have various "simulation" classes (either
    generated directly as here, or subclassed for more complication
    simulation objects as in OpenMM), but they all quack the same when it
    comes to things the DynamicsEngine calls on them for

    """

    base_snapshot_type = Snapshot

    default_options = {
        'n_frames_max': 5000,
        'nsteps_per_frame': 10
    }

    def __init__(self, msm, time_step=None, options=None):
        snapshot_dimensions = {
        }

        descriptor = SnapshotDescriptor.construct(
            snapshot_class=Snapshot,
            snapshot_dimensions=snapshot_dimensions
        )

        super(MSMEngine, self).__init__(
            options=options,
            descriptor=descriptor
        )

        self._state = None

        self.msm = msm
        if time_step is not None:
            self.time_step = time_step
        else:
            self.time_step = 1


        self._effective_msm = np.linalg.matrix_power(
            self.msm,
            self.nsteps_per_frame
        )

        self._n_states = len(msm)

    def to_dict(self):
        return {
            'options': self.options,
            'msm': self.msm,
            'time_step': self.time_step
        }

    @property
    def nsteps_per_frame(self):
        return self.options['nsteps_per_frame']

    @nsteps_per_frame.setter
    def nsteps_per_frame(self, value):
        self.options['nsteps_per_frame'] = value

    @property
    def snapshot_timestep(self):
        return self.nsteps_per_frame * self.time_step

    @property
    def state(self):
        return self._state

    @property
    def n_state(self):
        return self._n_states

    @property
    def current_snapshot(self):
        state = self._state
        return Snapshot(
            state=state,
            engine=self
        )

    @current_snapshot.setter
    def current_snapshot(self, snap):
        self.check_snapshot_type(snap)
        self._state = snap.state

    def generate_next_frame(self):
        self._state = np.random.choice(self._n_states, p=self.msm[self._state])
        return self.current_snapshot
