import logging

import numpy as np
from openpathsampling.engines import DynamicsEngine, SnapshotDescriptor

import hoomd

from .snapshot import Snapshot

logger = logging.getLogger(__name__)


class HOOMDEngine(DynamicsEngine):
    """HOOMD-blue dynamics engine based on ``hoomd`` system and integrator.

    The engine will create a :class:`hoomd.Simulation` instance
    and uses this to generate new frames.

    """

    _default_options = {
        "n_steps_per_frame": 10,
        "n_frames_max": 5000,
    }

    base_snapshot_type = Snapshot

    def __init__(self, simulation, options=None):
        """
        Parameters
        ----------
        sim : hoomd.Simulation
            A HOOMD Simulation object. This encompasses the system state,
            operations to apply when running, and the device to use for
            execution.
        options : dict
            A dictionary that provides additional settings for the OPS engine.
            Allowed keys are:

                'n_steps_per_frame' : int, default: 10
                    the number of integration steps per returned snapshot
                'n_frames_max' : int or None, default: 5000,
                    the maximal number of frames allowed for a returned
                    trajectory object

        Notes
        -----
        The `n_frames_max` does not limit Trajectory objects in length. It only
        limits the maximal length of returned trajectory objects when this
        engine is used.

        """

        self._simulation = simulation

        dimensions = {
            "n_atoms": self.sim.state.N_particles,
            "n_spatial": self.sim.state.box.dimensions,
        }

        descriptor = SnapshotDescriptor.construct(Snapshot, dimensions)

        super(HOOMDEngine, self).__init__(options=options, descriptor=descriptor)

        # Set no cached snapshot.
        self._current_snapshot = None
        self._current_momentum = None
        self._current_configuration = None
        self._current_box_vectors = None

        self._simulation = None
        self._n_dofs = None

    @property
    def simulation(self):
        """The `hoomd.Simulation` instance."""
        return self._simulation

    @property
    def mdtraj_topology(self):
        return self.topology.mdtraj

    @property
    def snapshot_timestep(self):
        integrator = self.simulation.operations.integrator
        if issubclass(integrator, hoomd.md.Integrator):
            dt = integrator.dt
        else:
            dt = 1
        return self.n_steps_per_frame * dt

    def _build_current_snapshot(self):
        hoomd_snapshot = self.simulation.state.get_snapshot()
        coordinates = np.array(hoomd_snapshot.particles.position)
        box_vectors = np.array(hoomd_snapshot.configuration.box.matrix)
        velocities = np.array(hoomd_snapshot.particles.velocity)

        return Snapshot.construct(
            coordinates=coordinates,
            box_vectors=box_vectors,
            velocities=velocities,
            engine=self,
        )

    @staticmethod
    def is_valid_snapshot(snapshot):
        return (not np.any(np.isnan(snapshot.coordinates._value))) and (
            not np.any(np.isnan(snapshot.velocities._value))
        )

    @property
    def current_snapshot(self):
        if self._current_snapshot is None:
            self._current_snapshot = self._build_current_snapshot()

        return self._current_snapshot

    def _changed(self):
        self._current_snapshot = None

    @current_snapshot.setter
    def current_snapshot(self, snapshot):
        self.check_snapshot_type(snapshot)

        if snapshot is not self._current_snapshot:
            hoomd_snapshot = self.simulation.state.get_snapshot(snapshot)
            if snapshot.coordinates is not None:
                hoomd_snapshot.particles.position = snapshot.coordinates

            if snapshot.box_vectors is not None:
                hoomd_snapshot.configuration.box = hoomd.Box.from_matrix(
                    snapshot.box_vectors
                )

            if snapshot.velocities is not None:
                hoomd_snapshot.particles.velocities = snapshot.velocities

            # After the updates cache the new snapshot
            if snapshot.engine is self:
                # no need for copy if this snap is from this engine
                self._current_snapshot = snapshot
            else:
                self._current_snapshot = self._build_current_snapshot()

    def generate_next_frame(self):
        self.simulation.run(self.n_steps_per_frame)
        self._current_snapshot = None
        return self.current_snapshot
