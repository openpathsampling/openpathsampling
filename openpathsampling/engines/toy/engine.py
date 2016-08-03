import numpy as np

from openpathsampling.engines import DynamicsEngine, SnapshotDescriptor
from snapshot import ToySnapshot as Snapshot


class ToyEngine(DynamicsEngine):
    """
    The trick is that we have various "simulation" classes (either
    generated directly as here, or subclassed for more complication
    simulation objects as in OpenMM), but they all quack the same when it
    comes to things the DynamicsEngine calls on them for

    """

    base_snapshot_type = Snapshot

    default_options = {
        'integ': None,
        'n_frames_max': 5000,
        'n_steps_per_frame': 10
    }

    def __init__(self, options, topology):
        if 'n_spatial' not in options:
            options['n_spatial'] = topology.n_spatial

        options['n_atoms'] = 1

        snapshot_dimensions = {
            'n_atoms': topology.n_atoms,
            'n_spatial': topology.n_spatial
        }

        descriptor = SnapshotDescriptor.construct(
            snapshot_class=Snapshot,
            snapshot_dimensions=snapshot_dimensions
        )

        super(ToyEngine, self).__init__(
            options=options,
            descriptor=descriptor
        )

        self.topology = topology

        self._mass = None
        self._minv = None

        self.positions = None
        self.velocities = None

        self._mass = np.array(topology.masses)
        self._pes = topology.pes
        self._minv = 1.0 / self._mass

    def to_dict(self):
        return {
            'options': self.options,
            'topology': self.topology
        }

    @property
    def pes(self):
        return self._pes

    @property
    def n_steps_per_frame(self):
        return self.options['n_steps_per_frame']

    @n_steps_per_frame.setter
    def n_steps_per_frame(self, value):
        self.options['n_steps_per_frame'] = value

    @property
    def mass(self):
        return self._mass

    @mass.setter
    def mass(self, value):
        self._mass = value
        self._minv = np.reciprocal(value)

    @property
    def snapshot_timestep(self):
        return self.n_steps_per_frame * self.integ.dt

    @property
    def current_snapshot(self):
        snap_pos = self.positions
        snap_vel = self.velocities
        return Snapshot(
            coordinates=np.array([snap_pos]),
            velocities=np.array([snap_vel]),
            engine=self
        )

    @current_snapshot.setter
    def current_snapshot(self, snap):
        self.check_snapshot_type(snap)

        coords = np.copy(snap.coordinates)
        vels = np.copy(snap.velocities)
        self.positions = coords[0]
        self.velocities = vels[0]

    def generate_next_frame(self):
        self.integ.step(self, self.n_steps_per_frame)
        return self.current_snapshot
