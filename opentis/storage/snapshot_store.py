from simtk.unit import Quantity, nanometers, kilojoules_per_mole, picoseconds
import numpy as np

from opentis.snapshot import Snapshot, Configuration, Momentum
from object_storage import ObjectStorage
from wrapper import savecache, loadcache
from opentis.trajectory import Trajectory

class SnapshotStorage(ObjectStorage):
    """
    An ObjectStorage for Snapshots. Allow to store Snapshots instances in a netcdf file.
    """

    def __init__(self, storage = None):
        super(SnapshotStorage, self).__init__(storage, Snapshot)

    @loadcache
    def load(self, idx=None):
        '''
        Load a snapshot from the storage.

        Parameters
        ----------

        Returns
        -------
        snapshot : Snapshot
            the loaded snapshot instance
        '''

        snapshot = Snapshot()

        configuration_idx = self.configuration_idx(idx)
        momentum_idx = self.momentum_idx(idx)
        momentum_reversed = self.momentum_reversed(idx)

        snapshot.configuration = self.storage.configuration.load(configuration_idx)
        snapshot.momentum = self.storage.momentum.load(momentum_idx)

        snapshot.reversed = momentum_reversed

        snapshot.idx[self.storage] = idx

        return snapshot

    def all(self):
        """
        Return a trajectory consisting of all (unordered) frames in the storage
        """
        t = Trajectory()
        count = self.count()
        for snapshot_idx in range(0,count):
            t.append(self.load(snapshot_idx))

        return t

    @savecache
    def save(self, snapshot, idx=None):
        """
        Add the current state of the snapshot in the database. If nothing has changed then the snapshot gets stored using the same snapshots as before. Saving lots of diskspace

        Parameters
        ----------
        snapshot : Snapshot()
            the snapshot to be saved
        idx : int or None
            if idx is not None the index will be used for saving in the storage. This might overwrite already existing trajectories!

        Notes
        -----
        This also saves all contained frames in the snapshot if not done yet.
        A single Snapshot object can only be saved once!
        """
        storage = self.storage

        if snapshot.configuration is not None:
            storage.configuration.save(snapshot.configuration)
            self.save_variable('snapshot_configuration_idx', idx, snapshot.configuration.idx[storage])
        else:
            self.save_variable('snapshot_configuration_idx', idx, -1)

        if snapshot.momentum is not None:
            storage.momentum.save(snapshot.momentum)
            self.save_variable('snapshot_momentum_idx', idx, snapshot.momentum.idx[storage])
        else:
            self.save_variable('snapshot_momentum_idx', idx, -1)

        self.save_variable('snapshot_momentum_reversed', idx, int(snapshot.reversed))


    def configuration_idx(self, idx):
        '''
        Load snapshot indices for snapshot with ID 'idx' from the storage

        ARGUMENTS

        idx (int) - ID of the snapshot

        Returns
        -------
        snapshot (list of int) - snapshot indices
        '''
        return int(self.load_variable('snapshot_configuration_idx', idx))

    def momentum_idx(self, idx):
        '''
        Load snapshot indices for snapshot with ID 'idx' from the storage

        ARGUMENTS

        idx (int) - ID of the snapshot

        Returns
        -------
        snapshot (list of int) - snapshot indices
        '''
        return int(self.load_variable('snapshot_momentum_idx', idx))


    def momentum_reversed(self, idx):
        '''
        Load snapshot with ID 'idx' from the storage and return a list of reversed indicators for the momenta

        Parameters
        ----------

        idx : int
            index of the snapshot

        Returns
        -------
        list of boolean
            list of boolean which frames in the snapshot are reversed
        '''
        return bool(self.load_variable('snapshot_momentum_reversed', idx))


    def _init(self):
        '''
        Initializes the associated storage to index configuration_indices in it
        '''
        super(SnapshotStorage, self)._init()

        self.init_variable('snapshot_configuration_idx', 'index', self.db,
                description="snapshot[snapshot] is the snapshot index (0..n_configuration-1) of snapshot 'snapshot'.")

        self.init_variable('snapshot_momentum_idx', 'index', self.db,
                description="snapshot[snapshot] is the snapshot index (0..n_momentum-1) 'frame' of snapshot 'snapshot'.")

        self.init_variable('snapshot_momentum_reversed', 'bool', self.db)

#=============================================================================================
# ORDERPARAMETER UTILITY FUNCTIONS
#=============================================================================================

    @property
    def op_configuration_idx(self):
        """
        Returns aa function that returns for an object of this storage the idx
        """
        def idx(obj):
            return obj.configuration.idx[self.storage]

        return idx

    @property
    def op_momentum_idx(self):
        """
        Returns aa function that returns for an object of this storage the idx
        """
        def idx(obj):
            return obj.momentum.idx[self.storage]

        return idx


class MomentumStorage(ObjectStorage):
    """
    An ObjectStorage for Momenta. Allows to store Momentum() instances in a netcdf file.
    """

    def __init__(self, storage = None):
        super(MomentumStorage, self).__init__(storage, Momentum)

    @savecache
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

        storage = self.storage

        # TODO: This should never be empty when it is called. Since a Momentum() instance has velocities
        # TODO: If it was load lazy then it should be registered as already saved and if a snapshot does not
        # TODO: have velocities then it does not have a Momentum object

        # Store momentum.
        if momentum._velocities is not None:
            storage.variables['momentum_velocities'][idx,:,:] = (momentum.velocities / (nanometers / picoseconds)).astype(np.float32)
        else:
            print 'ERROR : Momentum should not be empty'
        if momentum._kinetic_energy is not None:
            storage.variables['momentum_kinetic'][idx] = momentum.kinetic_energy / kilojoules_per_mole
        else:
            # TODO: No kinetic energy is not yet supported
            print 'Think about how to handle this. It should only be None if loaded lazy and in this case it will never be saved.'

        # Force sync to disk to avoid data loss.
        storage.sync()

    @loadcache
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

        if not (Momentum.load_lazy and lazy):
            v = storage.variables['momentum_velocities'][idx,:,:].astype(np.float32).copy()
            velocities = Quantity(v, nanometers / picoseconds)
            T = storage.variables['momentum_kinetic'][idx]
            kinetic_energy = Quantity(T, kilojoules_per_mole)

        else:
            velocities = None
            kinetic_energy = None

        momentum = Momentum(velocities=velocities, kinetic_energy=kinetic_energy)
        momentum.idx[storage] = idx

        return momentum

    def update_velocities(self, obj):
        storage = self.storage

        idx = obj.idx[self.storage]
        v = storage.variables['momentum_velocities'][idx,:,:].astype(np.float32).copy()
        velocities = Quantity(v, nanometers / picoseconds)

        obj.velocities = velocities

    def update_kinetic_energy(self, obj):
        storage = self.storage

        idx = obj.idx[self.storage]
        T = storage.variables['momentum_kinetic'][idx]
        kinetic_energy = Quantity(T, kilojoules_per_mole)

        obj.kinetic_energy = kinetic_energy


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

        super(MomentumStorage, self)._init()

        atoms = self.storage.atoms

        # define dimensions used in configuration_indices
        if 'atom' not in self.storage.dimensions:
            self.init_dimension('atom', atoms) # number of atoms in the simulated system

        if 'spatial' not in self.storage.dimensions:
            self.init_dimension('spatial', 3)  # number of spatial dimensions

        self.init_variable('momentum_velocities', 'float', (self.db, 'atom','spatial'), 'nm',
                description="velocities[momentum][atom][coordinate] are velocities of atom 'atom' in" +
                            " dimension 'coordinate' of momentum 'momentum'.")

        self.init_variable('momentum_kinetic', 'float', self.db)


    
class ConfigurationStorage(ObjectStorage):
    def __init__(self, storage = None):
        super(ConfigurationStorage, self).__init__(storage, Configuration)

    @savecache
    def save(self, configuration, idx = None):
        """
        Save positions, velocities, boxvectors and energies of current iteration to NetCDF file.

        Notes
        -----
        We need to allow for reversed configuration_indices to index memory. Would be nice
        """

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

    @loadcache
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

        if not (Configuration.load_lazy and lazy):
            x = storage.variables['configuration_coordinates'][idx,:,:].astype(np.float32).copy()
            coordinates = Quantity(x, nanometers)
            b = storage.variables['configuration_box_vectors'][idx]
            box_vectors = Quantity(b, nanometers)
            V = storage.variables['configuration_potential'][idx]
            potential_energy = Quantity(V, kilojoules_per_mole)
        else:
            coordinates = None
            box_vectors = None
            potential_energy = None

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

    def update_box_vectors(self, obj):
        storage = self.storage

        idx = obj.idx[self.storage]

        b = storage.variables['configuration_box_vectors'][idx]
        box_vectors = Quantity(b, nanometers)

        obj.box_vectors = box_vectors

    def update_potential_energy(self, obj):
        storage = self.storage

        idx = obj.idx[self.storage]

        V = storage.variables['configuration_potential'][idx]
        potential_energy = Quantity(V, kilojoules_per_mole)

        obj.potential_energy = potential_energy

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

    def snapshot_coordinates_as_array(self, idx, atom_indices=None):
        '''
        Returns a numpy array consisting of all coordinates of a snapshot

        Parameters
        ----------
        idx : int
            index of the snapshot to be loaded
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

        super(ConfigurationStorage, self)._init()

        atoms = self.storage.atoms

        # define dimensions used in configuration_indices
        if 'atom' not in self.storage.dimensions:
            self.init_dimension('atom', atoms) # number of atoms in the simulated system

        if 'spatial' not in self.storage.dimensions:
            self.init_dimension('spatial', 3)  # number of spatial dimensions

        self.init_variable('configuration_coordinates', 'float', (self.db, 'atom','spatial'), 'nm',
                description="coordinates[configuration][atom][coordinate] are coordinate of atom 'atom' " +
                            "in dimension 'coordinate' of configuration 'configuration'.")

        self.init_variable('configuration_box_vectors', 'float', (self.db, 'spatial'))

        self.init_variable('configuration_potential', 'float', self.db)
