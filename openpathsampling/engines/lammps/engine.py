import numpy as np
import simtk.unit as u

import openpathsampling as paths
from openpathsampling.engines import DynamicsEngine
from snapshot import Snapshot

from lammps import lammps


class LammpsEngine(DynamicsEngine):
    """OpenMM dynamics engine based on a openmmtools.testsystem object.

    This is only to allow to use the examples from openmmtools.testsystems
    within this framework
    """

    # units = {
    #     'length': u.dimensionless,
    #     'velocity': u.dimensionless,
    #     'energy': u.dimensionless
    # }

    _default_options = {
        'n_steps_per_frame': 10,
        'n_frames_max': 5000
    }

    def __init__(self, inputs, template=None, options=None):

        self.inputs = inputs

        # Create new lammps instance
        self._lmp = lammps()

        # Execute the give script
        commands = inputs.splitlines()

        for command in commands:
            self._lmp.command(command)

        # self.command('compute thermo_ke all ke')
        self.command('run 1')

        if template is None:
            template = self._get_snapshot(True)

        dimensions = {
            'n_atoms': template.coordinates.shape[0],
            'n_spatial': template.coordinates.shape[1]
        }

        descriptor = paths.engines.SnapshotDescriptor.construct(
            Snapshot,
            dimensions
        )

        super(LammpsEngine, self).__init__(
            options=options,
            descriptor=descriptor

        )

        # set no cached snapshot, means it will be constructed
        # from the openmm context
        self._current_snapshot = None
        self._current_momentum = None
        self._current_configuration = None
        self._current_box_vectors = None

        # TODO: so far we will always have an initialized system which we should change somehow
        self.initialized = True

    def command(self, *args, **kwargs):
        self._lmp.command(*args, **kwargs)

    def create(self):
        """
        Create the final OpenMMEngine

        """
        self.initialized = True

    def _get_snapshot(self, topology=None):
        lmp = self._lmp
        x = lmp.gather_atoms("x", 1, 3)
        v = lmp.gather_atoms("v", 1, 3)
        # pe = lmp.extract_compute('thermo_pe', 0, 0)
        # ke = lmp.extract_compute('thermo_ke', 0, 0)
        xlo = lmp.extract_global("boxxlo", 1)
        xhi = lmp.extract_global("boxxhi", 1)
        ylo = lmp.extract_global("boxylo", 1)
        yhi = lmp.extract_global("boxyhi", 1)
        zlo = lmp.extract_global("boxzlo", 1)
        zhi = lmp.extract_global("boxzhi", 1)
        xy = lmp.extract_global("xy", 1)
        xz = lmp.extract_global("xz", 1)
        yz = lmp.extract_global("yz", 1)
        bv = np.array(
            [[xhi - xlo, 0.0, 0.0], [xy, yhi - ylo, 0.0], [xz, yz, zhi - zlo]])
        n_atoms = lmp.get_natoms()
        # n_spatial = len(x) / n_atoms

        snapshot = Snapshot.construct(
            engine=self,
            coordinates=np.ctypeslib.array(x).reshape(
                (n_atoms, -1)) * u.nanometers,
            box_vectors=bv * u.nanometers,
            velocities=np.ctypeslib.array(v).reshape(
                (n_atoms, -1)) * u.nanometers / u.picoseconds
        )

        return snapshot

    def _put_coordinates(self, nparray):
        lmp = self._lmp
        lmparray = np.ctypeslib.as_ctypes(nparray.ravel())
        lmp.scatter_atoms("x", 1, 3, lmparray)

    def _put_velocities(self, nparray):
        lmp = self._lmp
        lmparray = np.ctypeslib.as_ctypes(nparray.ravel())
        lmp.scatter_atoms("v", 1, 3, lmparray)

    @property
    def lammps(self):
        return self._lmp

    def to_dict(self):
        return {
            'inputs': self.inputs,
            'template': self.template,
            'options': self.options
        }

    @classmethod
    def from_dict(cls, dct):
        inputs = dct['inputs']
        template = dct['template']
        options = dct['options']

        return LammpsEngine(
            inputs=inputs,
            template=template,
            options=options
        )

    @property
    def snapshot_timestep(self):
        return self.n_steps_per_frame * self.options['timestep']

    def _build_current_snapshot(self):
        return self._get_snapshot()

    @property
    def current_snapshot(self):
        if self._current_snapshot is None:
            self._current_snapshot = self._build_current_snapshot()

        return self._current_snapshot

    def _changed(self):
        self._current_snapshot = None

    @current_snapshot.setter
    def current_snapshot(self, snapshot):
        if snapshot is not self._current_snapshot:
            if snapshot.configuration is not None:
                if self._current_snapshot is None or snapshot.configuration is not self._current_snapshot.configuration:
                    # new snapshot has a different configuration so update
                    self._put_coordinates(snapshot.coordinates)

            if snapshot.momentum is not None:
                if self._current_snapshot is None or snapshot.momentum is not self._current_snapshot.momentum or snapshot.is_reversed != self._current_snapshot.is_reversed:
                    self._put_velocities(snapshot.velocities)

            # After the updates cache the new snapshot
            self._current_snapshot = snapshot

    def run(self, steps):
        self._lmp.command('run ' + str(steps))

    def generate_next_frame(self):
        self.run(self.n_steps_per_frame)
        self._current_snapshot = None
        return self.current_snapshot

    @property
    def momentum(self):
        return self.current_snapshot.momentum

    @property
    def configuration(self):
        return self.current_snapshot.configuration
