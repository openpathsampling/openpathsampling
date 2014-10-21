from snapshot import Snapshot, Configuration, Momentum
from object_storage import ObjectStorage
from simtk.unit import Quantity, nanometers, kilojoules_per_mole, picoseconds
import numpy as np

class SnapshotStorage(ObjectStorage):
    """
    An ObjectStorage for Snapshots. Allow to store Snapshots instances in a netcdf file.
    """

    def __init__(self, storage = None):
        super(SnapshotStorage, self).__init__(storage, Snapshot)

    def save(self, snapshot, idx_configuration = None, idx_momentum = None):
        """
        Save positions, velocities, boxvectors and energies of current iteration to NetCDF file.

        Parameters
        ----------
        snapshot : Snapshot()
            the actual snapshot object to be saved.
        idx_configuration : int or None
            if not None the configuration is saved using the index specified. Might overwrite an existing configuration.
        idx_momentum : int or None
            if not None the momentum is saved using the index specified. Might overwrite an existing momentum.
        """

        self.storage.configuration.save(snapshot.configuration, idx_configuration)
        self.storage.momentum.save(snapshot.momentum, idx_momentum)

    def load(self, idx_configuration = None, idx_momentum = None, reversed = False, lazy=None):
        '''
        Load a snapshot from the storage.

        Parameters
        ----------
        idx_configuration : int
            index of the configuration in the database 'idx' > 0
        idx:momentum : int
            index of the momentum in the database 'idx' > 0
        reversed : boolean
            if True the momenta are treated as reversed

        Returns
        -------
        snapshot : Snapshot
            the loaded snapshot instance
        '''

        #TODO: Check, for some reason some idx are given as numpy.in32 and netcdf4 is not compatible with indices given in this format!!!!!

        snapshot = Snapshot()
        snapshot.reversed = bool(reversed)

        if lazy is None:
            lazy_configuration = False
            lazy_momentum = True
        else:
            lazy_configuration = lazy
            lazy_momentum = lazy

        if idx_configuration is not None:
            idx_c = int(idx_configuration)
            snapshot.configuration = self.storage.configuration.load(idx_c, lazy=lazy_configuration)

        if idx_momentum is not None:
            idx_m = int(idx_momentum)
            snapshot.momentum = self.storage.momentum.load(idx_m, lazy=lazy_momentum)

        return snapshot

    def _init(self):
        pass

    def free(self):
        return 1

    def count(self):
        return 0

    def first(self):
        return None

    def last(self):
        return None

    def get(self, indices):
        return None

class MomentumStorage(ObjectStorage):
    """
    An ObjectStorage for Momenta. Allows to store Momentum() instances in a netcdf file.
    """

    def __init__(self, storage = None):
        super(MomentumStorage, self).__init__(storage, Momentum)

    def save(self, momentum, idx = None):
        """
        Save velocities and kinetic energies of current iteration to NetCDF file.

        Parameters
        ----------
        momentum : Momentum()
            the actual Momentum() instance to be saved.
        idx : int or None
            if not None `idx`is used as the index to index the Momentum() instance. Might overwrite existing Momentum in the database.
        """

        idx = super(MomentumStorage, self).index(momentum, idx)

        if idx is not None:
            storage = self.storage

            # Store momentum.
            storage.variables['momentum_velocities'][idx,:,:] = (momentum.velocities / (nanometers / picoseconds)).astype(np.float32)
            if momentum.kinetic_energy is not None:
                storage.variables['momentum_kinetic'][idx] = momentum.kinetic_energy / kilojoules_per_mole

            # Force sync to disk to avoid data loss.
            storage.sync()

        return


    def load(self, idx, lazy=True):
        '''
        Load a momentum from the storage

        Parameters
        ----------
        idx : int
            index of the momentum in the database 'idx' > 0

        Returns
        -------
        Momentum()
            the loaded momentum instance
        '''

        storage = self.storage

        #TODO: Check, for some reason some idx are given as numpy.in32 and netcdf4 is not compatible with indices given in this format!!!!!
        idx = int(idx)

        if not (Momentum.load_lazy and lazy):
            v = storage.variables['momentum_velocities'][idx,:,:].astype(np.float32).copy()
            velocities = Quantity(v, nanometers / picoseconds)
        else:
            velocities = None

        T = storage.variables['momentum_kinetic'][idx]
        kinetic_energy = Quantity(T, kilojoules_per_mole)

        momentum = Momentum(velocities=velocities, kinetic_energy=kinetic_energy)
        momentum.idx[storage] = idx

        return momentum

    def update_velocities(self, obj):
        storage = self.storage

        #TODO: Check, for some reason some idx are given as numpy.in32 and netcdf4 is not compatible with indices given in this format!!!!!
        idx = obj.idx[self.storage]

        v = storage.variables['momentum_velocities'][idx,:,:].astype(np.float32).copy()
        velocities = Quantity(v, nanometers / picoseconds)

        obj.velocities = velocities

    def velocities_as_numpy(self, frame_indices=None, atom_indices=None):
        """
        Return a block of stored velocities in the database as a numpy array.

        Parameters
        ----------
        frame_indices : list of int or None
            the indices of Momentum objects to be retrieved from the database. If `None` is specified then all indices are returned!
        atom_indices : list of int of None
            if not None only the specified atom_indices are returned. Might speed up reading a lot.
        """

        if frame_indices is None:
            frame_indices = slice(None)

        if atom_indices is None:
            atom_indices = slice(None)

        return self.variables['momentum_velocities'][frame_indices,atom_indices,:].astype(np.float32).copy()

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

    def _init(self):
        '''
        Initializes the associated storage to index momentums in it
        '''

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
        
    
class ConfigurationStorage(ObjectStorage):
    def __init__(self, storage = None):
        super(ConfigurationStorage, self).__init__(storage, Configuration)

    def save(self, configuration, idx = None):
        """
        Save positions, velocities, boxvectors and energies of current iteration to NetCDF file.

        Notes
        -----
        We need to allow for reversed configuration_indices to index memory. Would be nice
        """

        idx = super(ConfigurationStorage, self).index(configuration, idx)
        if idx is not None:
            storage = self.storage

            # Store configuration.
            storage.variables['configuration_coordinates'][idx,:,:] = (configuration.coordinates / nanometers).astype(np.float32)

            if configuration.potential_energy is not None:
                storage.variables['configuration_potential'][idx] = configuration.potential_energy / kilojoules_per_mole
    #            storage.variables['configuration_box_vectors'][idx,:] = (self.box_vectors / nanometers).astype(np.float32)

            # Force sync to disk to avoid data loss.
            storage.sync()


    def coordinates_as_numpy(self, frame_indices=None, atom_indices=None, storage = None):
        return self.coordinates_as_numpy(self, frame_indices, atom_indices)

    def get(self, indices):
        return [ self.load(idx) for idx in indices ]

    def load(self, idx, lazy=False):
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
        if not (Configuration.load_lazy and lazy):
            x = storage.variables['configuration_coordinates'][idx,:,:].astype(np.float32).copy()
            coordinates = Quantity(x, nanometers)
        else:
            coordinates = None

        b = storage.variables['configuration_box_vectors'][idx]
        box_vectors = Quantity(b, nanometers)
        V = storage.variables['configuration_potential'][idx]
        potential_energy = Quantity(V, kilojoules_per_mole)

        configuration = Configuration(coordinates=coordinates, box_vectors = box_vectors, potential_energy=potential_energy)
        configuration.idx[storage] = idx

        configuration.topology = self.storage.topology

        return configuration

    def update_coordinates(self, obj):
        storage = self.storage

        idx = obj.idx[self.storage]

        x = storage.variables['configuration_coordinates'][idx,:,:].astype(np.float32).copy()
        coordinates = Quantity(x, nanometers)

        obj.coordinates = coordinates

    def coordinates_as_numpy(self, frame_indices=None, atom_indices=None):

        if frame_indices is None:
            frame_indices = slice(None)

        if atom_indices is None:
            atom_indices = slice(None)

        return self.storage.variables['configuration_coordinates'][frame_indices,atom_indices,:].astype(np.float32).copy()

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
        Initializes the associated storage to index configuration_indices in it
        '''
        # index associated storage in class variable for all configuration instances to access

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
