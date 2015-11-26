from openpathsampling.netcdfplus import ObjectStore
from openpathsampling.tools import units_from_snapshot

from openpathsampling.snapshot import Configuration

import numpy as np

_variables = ['configuration']


def _init(store):
    store.storage.create_store('configurations', ConfigurationStore())

    store.create_variable('configuration', 'lazyobj.configurations',
                        description="the snapshot index (0..n_configuration-1) of snapshot '{idx}'.",
                        chunksizes=(1,)
                        )


class ConfigurationStore(ObjectStore):
    """
    An ObjectStore for Configuration. Allows to store Configuration() instances in a netcdf file.
    """
    def __init__(self):
        super(ConfigurationStore, self).__init__(Configuration, json=False)

    def _save(self, configuration, idx):
        # Store configuration.
        self.vars['coordinates'][idx] = configuration.coordinates

        if configuration.potential_energy is not None:
            self.vars['potential_energy'][idx] = configuration.potential_energy

        if configuration.box_vectors is not None:
            self.vars['box_vectors'][idx] = configuration.box_vectors

    def get(self, indices):
        return [self.load(idx) for idx in indices]

    def _load(self, idx):
        coordinates = self.vars["coordinates"][idx]
        box_vectors = self.vars["box_vectors"][idx]
        potential_energy = self.vars["potential_energy"][idx]

        configuration = Configuration(coordinates=coordinates, box_vectors=box_vectors,
                                      potential_energy=potential_energy)
        configuration.topology = self.storage.topology

        return configuration

    def coordinates_as_numpy(self, frame_indices=None, atom_indices=None):
        """
        Return the atom coordinates in the storage for given frame indices
        and atoms

        Parameters
        ----------
        frame_indices : list of int or None
            the frame indices to be included. If None all frames are returned
        atom_indices : list of int or None
            the atom indices to be included. If None all atoms are returned

        Returns
        -------
        numpy.array, shape=(n_frames, n_atoms)
            the array of atom coordinates in a float32 numpy array

        """
        if frame_indices is None:
            frame_indices = slice(None)

        if atom_indices is None:
            atom_indices = slice(None)

        return self.storage.variables[self.prefix + '_coordinates'][frame_indices, atom_indices, :].astype(
            np.float32).copy()

    def _init(self):
        super(ConfigurationStore, self)._init()
        n_atoms = self.storage.n_atoms
        n_spatial = self.storage.n_spatial

        units = units_from_snapshot(self.storage._template)

        self.create_variable('coordinates', 'numpy.float32',
                           dimensions=('atom', 'spatial'),
                           description="coordinate of atom '{ix[1]}' in dimension " +
                                       "'{ix[2]}' of configuration '{ix[0]}'.",
                           chunksizes=(1, n_atoms, n_spatial),
                           simtk_unit=units['length']
                           )

        self.create_variable('box_vectors', 'numpy.float32',
                           dimensions=('spatial', 'spatial'),
                           chunksizes=(1, n_spatial, n_spatial),
                           simtk_unit=units['length']
                           )

        self.create_variable('potential_energy', 'float',
                           chunksizes=(1,),
                           simtk_unit=units['energy']
                           )

