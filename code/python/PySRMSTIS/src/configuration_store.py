from simtk.unit import Quantity, nanometers, kilojoules_per_mole
import numpy as np

from snapshot import Configuration

from storage_utils import setstorage

from functools import wraps
from storage_utils import ObjectStorage

class ConfigurationStorage(ObjectStorage):
    def __init__(self, storage = None):
        super(ConfigurationStorage, self).__init__(storage, Configuration)

    def save(self, configuration, idx = None):
        """
        Save positions, velocities, boxvectors and energies of current iteration to NetCDF file.

        Notes
        -----
        We need to allow for reversed configuration_indices to save memory. Would be nice
        """

        idx = super(ConfigurationStorage, self).save(configuration, idx)
        if idx is not None:
            storage = self.storage

            # Store configuration.
            storage.variables['configuration_coordinates'][idx,:,:] = (configuration.coordinates / nanometers).astype(np.float32)

            if configuration.potential_energy is not None:
                storage.variables['configuration_potential'][idx] = configuration.potential_energy / kilojoules_per_mole
    #            storage.variables['configuration_box_vectors'][idx,:] = (self.box_vectors / nanometers).astype(np.float32)

            # Force sync to disk to avoid data loss.
            storage.sync()


    @wraps(setstorage)
    def coordinates_as_numpy(self, frame_indices=None, atom_indices=None, storage = None):
        return self.coordinates_as_numpy(self, frame_indices, atom_indices)

    @wraps(setstorage)
    def get(self, indices):
        return [ self.load(idx) for idx in indices ]

    @wraps(setstorage)
    def load(self, idx):
        '''
        Load a configuration from the storage

        Parameters
        ----------
        idx : int
            index of the configuration in the database 'idx' > 0

        Returns
        -------
        configuration : configuration
            the configuration
        '''

        storage = self.storage

        #TODO: Check, for some reason some idx are given as numpy.in32 and netcdf4 is not compatible with indices given in this format!!!!!
        idx = int(idx)

        #TODO: Use newest simtk.units since there was an inconcistance with the new numpy
        x = storage.variables['configuration_coordinates'][idx,:,:].astype(np.float32).copy()
        coordinates = Quantity(x, nanometers)
        b = storage.variables['configuration_box_vectors'][idx]
        box_vectors = Quantity(b, nanometers)
        V = storage.variables['configuration_potential'][idx]
        potential_energy = Quantity(V, kilojoules_per_mole)

        configuration = Configuration(coordinates=coordinates, box_vectors = box_vectors, potential_energy=potential_energy)
        configuration.idx[storage] = idx

        configuration.topology = self.storage.topology

        return configuration

    @wraps(setstorage)
    def coordinates_as_numpy(self, frame_indices=None, atom_indices=None):

        if frame_indices is None:
            frame_indices = slice(None)

        if atom_indices is None:
            atom_indices = slice(None)

        return self.storage.variables['configuration_coordinates'][frame_indices,atom_indices,:].astype(np.float32).copy()

    @wraps(setstorage)
    def coordinates_as_array(self, frame_indices=None, atom_indices=None):
        '''
        Returns a numpy array consisting of all coordinates at the given indices

        Parameters
        ----------
        frame_indices : list of int
            configuration indices to be loaded
        atom_indices : list of int
            selects only the atoms to be returned. If None (Default) all atoms will be selected

        Returns
        -------
        numpy.ndarray, shape = (l,n)
            returns an array with `l` the number of frames and `n` the number of atoms
        '''

        return self.coordinates_as_numpy(frame_indices, atom_indices)

    def trajectory_coordinates_as_array(self, idx, atom_indices=None):
        '''
        Returns a numpy array consisting of all coordinates of a trajectory

        Parameters
        ----------
        idx : int
            index of the trajectory to be loaded
        atom_indices : list of int
            selects only the atoms to be returned. If None (Default) all atoms will be selected


        Returns
        -------
        numpy.ndarray, shape = (l,n)
            returns an array with `l` the number of frames and `n` the number of atoms
        '''

        frame_indices = self.configuration_indices(idx)
        return self.coordinates_as_array(frame_indices, atom_indices)

    def _init(self):
        '''
        Initializes the associated storage to save configuration_indices in it
        '''
        # save associated storage in class variable for all configuration instances to access

#        ncgrp = storage.createGroup('configuration')

        super(ConfigurationStorage, self)._init()

        ncgrp = self.storage
        atoms = self.storage.atoms

        # define dimensions used in configuration_indices
        if 'atom' not in ncgrp.dimensions:
            ncgrp.createDimension('atom', atoms)    # number of atoms in the simulated system

        if 'spatial' not in ncgrp.dimensions:
            ncgrp.createDimension('spatial', 3) # number of spatial dimensions

#        print 'NAME', self.idx_dimension

        # define variables for configuration_indices
        ncvar_configuration_coordinates          = ncgrp.createVariable('configuration_coordinates', 'f', (self.idx_dimension,'atom','spatial'))
        ncvar_configuration_box_vectors          = ncgrp.createVariable('configuration_box_vectors', 'f', (self.idx_dimension, 'spatial'))
        ncvar_configuration_potential            = ncgrp.createVariable('configuration_potential',   'f', self.idx_dimension)

        # Define units for configuration variables.
        setattr(ncvar_configuration_coordinates, 'units', 'nm')
        setattr(ncvar_configuration_box_vectors,  'units', 'nm')
        setattr(ncvar_configuration_potential,   'units', 'kJ/mol')

        # Define long (human-readable) names for variables.
        setattr(ncvar_configuration_coordinates,   "long_name", "coordinates[configuration][atom][coordinate] are coordinate of atom 'atom' in dimension 'coordinate' of configuration 'configuration'.")
