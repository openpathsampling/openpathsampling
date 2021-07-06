import numpy as np

from openpathsampling.engines import DynamicsEngine, SnapshotDescriptor
from .snapshot import ToySnapshot as Snapshot


class ToyEngine(DynamicsEngine):
    """Engine for toy models. Mostly used for 2D examples.

    Parameters
    ----------
    options : dict
        A dictionary providing additional settings. Keys can be

            'integ' : :class:`.ToyIntegrator`
                the integrator for this engine
            'n_frames_max' : int
                the maximum number of frames allowed for a returned
                trajectory, default is 5000
            'n_steps_per_frame' : int
                number of integration steps per returned snapshot, default
                is 10.

    topology : :class:`.ToyTopology`
        object which includes masses, potential energy surface, and the
        dimensions n_atoms and n_spatial; plays a role similar to a topology
        in molecular mechanics

    Attributes
    ----------
    pes : :class:`.PES`
        potential energy surface
    mass : array-like
        mass of each atom
    snapshot_timestep : float
        time step between reported snapshots
    current_snapshot : :class:`.Snapshot`
        the current state of the system, as a snapshot
    """

    base_snapshot_type = Snapshot
    ignore_linear_momentum = True

    _default_options = {
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
        for i in range(self.n_steps_per_frame):
            self.integ.step(sys=self)
        return self.current_snapshot

    def n_degrees_of_freedom(self):
        topol = self.topology
        return topol.n_atoms * topol.n_spatial

    def has_constraints(self):
        return False
