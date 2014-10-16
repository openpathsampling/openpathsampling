__author__ = 'jan-hendrikprinz'

from simtk.unit import Quantity, nanometers, kilojoules_per_mole, picoseconds

import numpy as np

from snapshot import Momentum

from storage_utils import setstorage

from functools import wraps
from storage_utils import ObjectStorage

class MomentumStorage(ObjectStorage):
    def __init__(self, storage = None):
        super(MomentumStorage, self).__init__(storage, Momentum)

    @wraps(setstorage)
    def save(self, momentum, idx = None):
        """
        Save velocities and kinetic energies of current iteration to NetCDF file.
        """

        idx = super(MomentumStorage, self).save(momentum, idx)

        if idx is not None:
            storage = self.storage

            # Store momentum.
            storage.variables['momentum_velocities'][idx,:,:] = (momentum.velocities / (nanometers / picoseconds)).astype(np.float32)
            if momentum.kinetic_energy is not None:
                storage.variables['momentum_kinetic'][idx] = momentum.kinetic_energy / kilojoules_per_mole

            # Force sync to disk to avoid data loss.
            storage.sync()

        return


    def load(self, idx):
        '''
        Load a momentum from the storage

        Parameters
        ----------
        idx : int
            index of the momentum in the database 'idx' > 0

        Returns
        -------
        momentum : self
            the momentum
        '''

        storage = self.storage

        #TODO: Check, for some reason some idx are given as numpy.in32 and netcdf4 is not compatible with indices given in this format!!!!!
        idx = int(idx)

        v = storage.variables['momentum_velocities'][idx,:,:].astype(np.float32).copy()
        velocities = Quantity(v, nanometers / picoseconds)
        T = storage.variables['momentum_kinetic'][idx]
        kinetic_energy = Quantity(T, kilojoules_per_mole)

        momentum = Momentum(velocities=velocities, kinetic_energy=kinetic_energy)
        momentum.idx[storage] = idx

        return momentum

    @wraps(setstorage)
    def velocities_as_numpy(self, frame_indices=None, atom_indices=None):

        if frame_indices is None:
            frame_indices = slice(None)

        if atom_indices is None:
            atom_indices = slice(None)

        return self.variables['momentum_velocities'][frame_indices,atom_indices,:].astype(np.float32).copy()

    @wraps(setstorage)
    def velocities_as_array(self, frame_indices=None, atom_indices=None):
        '''
        Returns a numpy array consisting of all velocities at the given indices

        Parameters
        ----------
        frame_indices : list of int
            momenta indices to be loaded
        atom_indices : list of int
            selects only the atoms to be returned. If None (Default) all atoms will be selected


        Returns
        -------
        numpy.ndarray, shape = (l,n)
            returns an array with `l` the number of frames and `n` the number of atoms
        '''

        return self.velocities_as_numpy(frame_indices, atom_indices)

    @wraps(setstorage)
    def trajectory_velocities_as_array(self, idx, atom_indices=None):
        '''
        Returns a numpy array consisting of all velocities of a trajectory

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

        frame_indices = self.load_configuration_indices(idx)
        return self.coordinates_as_array(frame_indices, atom_indices)

    @wraps(setstorage)
    def _init(self):
        '''
        Initializes the associated storage to save momentums in it
        '''
        # save associated storage in class variable for all self instances to access

#        ncgrp = storage.createGroup('momentum')

        ncgrp = self.storage

        atoms = self.storage.atoms

        # define dimensions used in momentums
        ncgrp.createDimension('momentum', 0)                       # unlimited number of momentums
        if 'atom' not in ncgrp.dimensions:
            ncgrp.createDimension('atom', atoms)    # number of atoms in the simulated system

        if 'spatial' not in ncgrp.dimensions:
            ncgrp.createDimension('spatial', 3) # number of spatial dimensions

        # define variables for momentums
        ncvar_momentum_velocities           = ncgrp.createVariable('momentum_velocities',  'f', ('momentum','atom','spatial'))
        ncvar_momentum_kinetic              = ncgrp.createVariable('momentum_kinetic',     'f', ('momentum'))

        # Define units for momentum variables.
        setattr(ncvar_momentum_velocities,  'units', 'nm/ps')
        setattr(ncvar_momentum_kinetic,     'units', 'kJ/mol')

        # Define long (human-readable) names for variables.
        setattr(ncvar_momentum_velocities,    "long_name", "velocities[momentum][atom][coordinate] are velocities of atom 'atom' in dimension 'coordinate' of momentum 'momentum'.")