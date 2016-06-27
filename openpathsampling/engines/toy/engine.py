import numpy as np

from openpathsampling.engines import DynamicsEngine
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
                      'integ' : None,
                      'n_frames_max' : 5000,
                      'nsteps_per_frame' : 10
    }

    def __init__(self, options, topology):
        if 'n_spatial' not in options:
            options['n_spatial'] = topology.n_spatial

        options['n_atoms'] = 1

        super(ToyEngine, self).__init__(
            options=options,
            topology=topology
        )

        self.mass = topology.masses
        self._pes = topology.pes

    @property
    def pes(self):
        return self._pes

    @property
    def nsteps_per_frame(self):
        return self.options['nsteps_per_frame']

    @nsteps_per_frame.setter
    def nsteps_per_frame(self, value):
        self.options['nsteps_per_frame'] = value

    @property
    def mass(self):
        return self._mass

    @mass.setter
    def mass(self, value):
        self._mass = value
        self._minv = np.reciprocal(value)

    @property
    def snapshot_timestep(self):
        return self.nsteps_per_frame * self.integ.dt

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
        self.integ.step(self, self.nsteps_per_frame)
        return self.current_snapshot
