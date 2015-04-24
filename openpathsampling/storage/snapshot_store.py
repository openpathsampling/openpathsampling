import numpy as np
from openpathsampling.snapshot import Snapshot, Configuration, Momentum
from object_storage import ObjectStore
from openpathsampling.trajectory import Trajectory
import simtk.unit as u

class SnapshotStore(ObjectStore):
    """
    An ObjectStore for Snapshots in netCDF files.
    """

    def __init__(self, storage = None):
        super(SnapshotStore, self).__init__(storage, Snapshot, json=False)

    def load(self, idx=None):
        '''
        Load a snapshot from the storage.

        Parameters
        ----------
        idx : int
            the integer index of the snapshot to be loaded

        Returns
        -------
        snapshot : Snapshot
            the loaded snapshot instance
        '''

        snapshot = Snapshot()

        configuration_idx = self.configuration_idx(idx)
        momentum_idx = self.momentum_idx(idx)
        momentum_reversed = self.momentum_reversed(idx)

        snapshot.configuration = self.storage.configurations.load(configuration_idx)
        snapshot.momentum = self.storage.momentum.load(momentum_idx)

        snapshot.reversed = momentum_reversed

        return snapshot

    def all(self):
        """
        Return a trajectory consisting of all (unordered) frames in the storage.

        Notes
        -----
        If you are interested in collectivevariables this is faster since it does not
        load the snapshots. Otherwise storage.snapshots is fine to get an
        iterator. Both should should be about the same speed.
        """
        #TODO: Might think about replacing the iterator with this since it is
        # faster for collectivevariables
        return Trajectory([ (self, idx) for idx in range(len(self)) ])

    def save(self, snapshot, idx=None):
        """
        Add the current state of the snapshot in the database.

        Parameters
        ----------
        snapshot : Snapshot()
            the snapshot to be saved
        idx : int or None
            if idx is not None the index will be used for saving in the storage.
            This might overwrite already existing trajectories!

        Notes
        -----
        This also saves all contained frames in the snapshot if not done yet.
        A single Snapshot object can only be saved once!
        """
        storage = self.storage

        if snapshot.configuration is not None:
            storage.configurations.save(snapshot.configuration)
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
        Load snapshot index for snapshot with ID 'idx' from the storage

        Parameters
        ----------
        idx : int
            index of the snapshot

        Returns
        -------
        list of int
            configuration indices
        '''
        return int(self.load_variable('snapshot_configuration_idx', idx))

    def momentum_idx(self, idx):
        '''
        Load momentum index for snapshot with ID 'idx' from the storage

        Parameters
        ----------
        idx : int
            index of the snapshot

        Returns
        -------
        int
            momentum indices
        '''
        return int(self.load_variable('snapshot_momentum_idx', idx))


    def momentum_reversed(self, idx):
        '''
        Load reversed boolean for snapshot with ID 'idx' from the storage

        Parameters
        ----------
        idx : int
            index of the snapshot

        Returns
        -------
        boolean
            boolean if the momentum of the snapshot is reversed
        '''
        return bool(self.load_variable('snapshot_momentum_reversed', idx))


    def _init(self):
        '''
        Initializes the associated storage to index configuration_indices in it
        '''
        super(SnapshotStore, self)._init()

        self.init_variable('snapshot_configuration_idx', 'index', self.db,
                description="snapshot[snapshot] is the snapshot index (0..n_configuration-1) of snapshot 'snapshot'.",
                chunksizes=(1, )
        )

        self.init_variable('snapshot_momentum_idx', 'index', self.db,
                description="snapshot[snapshot] is the snapshot index (0..n_momentum-1) 'frame' of snapshot 'snapshot'.",
                chunksizes=(1, )
                )

        self.init_variable('snapshot_momentum_reversed', 'bool', self.db, chunksizes=(1, ))

#=============================================================================================
# COLLECTIVE VARIABLE UTILITY FUNCTIONS
#=============================================================================================

    @property
    def op_configuration_idx(self):
        """
        Returns aa function that returns for an object of this storage the idx

        Returns
        -------
        function
            the function that returns the idx of the configuration
        """
        def idx(obj):
            return obj.configuration.idx[self.storage]

        return idx

    @property
    def op_momentum_idx(self):
        """
        Returns aa function that returns for an object of this storage the idx

        Returns
        -------
        function
            the function that returns the idx of the configuration

        """
        def idx(obj):
            return obj.momentum.idx[self.storage]

        return idx


class MomentumStore(ObjectStore):
    """
    An ObjectStore for Momenta. Allows to store Momentum() instances in a netcdf file.
    """

    def __init__(self, storage = None):
        super(MomentumStore, self).__init__(storage, Momentum, json=False, load_partial=True)

        # attach delayed loaders
        self.set_variable_partial_loading('velocities', self.update_velocities)
        self.set_variable_partial_loading('kinetic_energy', self.update_kinetic_energy)


    def save(self, momentum, idx = None):
        """
        Save velocities and kinetic energies.

        Parameters
        ----------
        momentum : Momentum()
            the actual Momentum() instance to be saved.
        idx : int or None
            if not None `idx`is used as the index to index the Momentum()
            instance. Might overwrite existing Momentum in the database.
        """

        storage = self.storage

        if momentum.velocities is not None:
            if hasattr(momentum.velocities, 'unit'):
                storage.variables['momentum_velocities'][idx,:,:] = (momentum.velocities / self.storage.units["momentum_velocities"]).astype(np.float32)
            else:
                # in this case we blindly assume that the units are correct
                storage.variables['momentum_velocities'][idx,:,:] = momentum.velocities.astype(np.float32)
        else:
            print 'ERROR : Momentum should not be empty'

        if momentum.kinetic_energy is not None:
            storage.variables['momentum_kinetic'][idx] = momentum.kinetic_energy / self.storage.units["momentum_kinetic"]
        else:
            # TODO: No kinetic energy is not yet supported
            print 'Think about how to handle this. It should only be None if loaded lazy and in this case it will never be saved.'

        # Force sync to disk to avoid data loss.
        # storage.sync()

    def load(self, idx):
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

        v = storage.variables['momentum_velocities'][idx,:,:].astype(np.float32).copy()
        velocities = u.Quantity(v, self.storage.units["momentum_velocities"])
        T = storage.variables['momentum_kinetic'][idx]
        kinetic_energy = u.Quantity(T, self.storage.units["momentum_kinetic"])
        momentum = Momentum(velocities=velocities, kinetic_energy=kinetic_energy)

        return momentum

    def load_empty(self, idx):
        momentum = Momentum()
        del momentum.velocities
        del momentum.kinetic_energy
        return momentum

    @staticmethod
    def update_velocities(obj):
        """
        Update/Load the velocities in the given obj from the attached storage

        Parameters
        ----------
        obj : Momentum
            The Momentum object to be updated

        """
        storage = obj._origin

        idx = obj.idx[storage]
        v = storage.variables['momentum_velocities'][idx,:,:].astype(np.float32).copy()
        velocities = u.Quantity(v, storage.units["momentum_velocities"])

        obj.velocities = velocities

    @staticmethod
    def update_kinetic_energy(obj):
        """
        Update/Load the kinetic_energy in the given obj from the attached storage

        Parameters
        ----------
        obj : Momentum
            The Momentum object to be updated

        """
        storage = obj._origin

        idx = obj.idx[storage]
        T = storage.variables['momentum_kinetic'][idx]
        kinetic_energy = u.Quantity(T, storage.units["momentum_kinetic"])

        obj.kinetic_energy = kinetic_energy


    def velocities_as_numpy(self, frame_indices=None, atom_indices=None):
        """
        Return a block of stored velocities in the database as a numpy array.

        Parameters
        ----------
        frame_indices : list of int or None
            the indices of Momentum objects to be retrieved from the database.
            If `None` is specified then all indices are returned!
        atom_indices : list of int of None
            if not None only the specified atom_indices are returned. Might
            speed up reading a lot.
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
            selects only the atoms to be returned. If None (Default) all atoms
            will be selected


        Returns
        -------
        numpy.ndarray, shape = (l,n)
            returns an array with `l` the number of frames and `n` the number
            of atoms
        '''

        return self.velocities_as_numpy(frame_indices, atom_indices)

    def _init(self):
        '''
        Initializes the associated storage to index momentums in it
        '''

        super(MomentumStore, self)._init()

        n_atoms = self.storage.n_atoms
        n_spatial = self.storage.n_spatial

        self.init_variable('momentum_velocities', 'float',
                (self.db, 'atom','spatial'),
                self.dimension_units['velocity'],
                description="velocities[momentum][atom][coordinate] are " +
                            "velocities of atom 'atom' in dimension " +
                            "'coordinate' of momentum 'momentum'.",
                chunksizes=(1,n_atoms,n_spatial))

        self.init_variable('momentum_kinetic', 'float', self.db,
                self.dimension_units['energy'],
                chunksizes=(1, ))


    
class ConfigurationStore(ObjectStore):
    def __init__(self, storage = None):
        super(ConfigurationStore, self).__init__(storage, Configuration, json=False, load_partial=True)

        # attach delayed loaders
        self.set_variable_partial_loading('coordinates', self.update_coordinates)
        self.set_variable_partial_loading('box_vectors', self.update_box_vectors)
        self.set_variable_partial_loading('potential_energy', self.update_potential_energy)

    def save(self, configuration, idx = None):
        storage = self.storage

        # Store configuration.
        storage.variables['configuration_coordinates'][idx,:,:] = (configuration.coordinates / self.storage.units["configuration_coordinates"]).astype(np.float32)

        if configuration.potential_energy is not None:
            storage.variables['configuration_potential'][idx] = configuration.potential_energy / self.storage.units["configuration_potential"]

        if configuration.box_vectors is not None:
            storage.variables['configuration_box_vectors'][idx,:,:] = (configuration.box_vectors / self.storage.units["configuration_box_vectors"]).astype(np.float32)

        # TODO: Add simple test if topologies match
        # if configuration.topology is not storage.topology:
        # log that topologies were different

        # Force sync to disk to avoid data loss.
        # storage.sync()


    def get(self, indices):
        return [ self.load(idx) for idx in indices ]

    def load(self, idx):
        storage = self.storage
        x = storage.variables['configuration_coordinates'][idx,:,:].astype(np.float32).copy()
        coordinates = u.Quantity(x, self.storage.units["configuration_coordinates"])
        b = storage.variables['configuration_box_vectors'][idx]
        box_vectors = u.Quantity(b, self.storage.units["configuration_box_vectors"])
        V = storage.variables['configuration_potential'][idx]
        potential_energy = u.Quantity(V, self.storage.units["configuration_potential"])

        configuration = Configuration(coordinates=coordinates, box_vectors = box_vectors, potential_energy=potential_energy)
        configuration.topology = self.storage.topology

#        print 'loaded normally'

        return configuration

    def load_empty(self, idx):
        """
        Loading function for partial loading. Constructs an empty Configuration
        object.

        Parameters
        ----------
        idx : int
            the integer index of the configuration to be loaded

        Returns
        -------
        Configuration
            an empty configuration object
        """
        configuration = Configuration()
        configuration.topology = self.storage.topology

        # if these still exist they will not be loaded using __getattr__
        del configuration.coordinates
        del configuration.box_vectors
        del configuration.potential_energy

        return configuration

    @staticmethod
    def update_coordinates(obj):
        """
        Update/Load the coordinates in the given obj from the attached storage

        Parameters
        ----------
        obj : Configuration
            the Configuration object to be updated

        """
        storage = obj._origin
        idx = obj.idx[storage]

        x = storage.variables['configuration_coordinates'][idx,:,:].astype(np.float32).copy()
        coordinates = u.Quantity(x, storage.units["configuration_coordinates"])

        obj.coordinates = coordinates

    @staticmethod
    def update_box_vectors(obj):
        """
        Update/Load the box_vectors in the given obj from the attached storage

        Parameters
        ----------
        obj : Configuration
            the Configuration object to be updated

        """
        storage = obj._origin
        idx = obj.idx[storage]

        b = storage.variables['configuration_box_vectors'][idx]
        box_vectors = u.Quantity(b, storage.units["configuration_box_vectors"])

        obj.box_vectors = box_vectors

    @staticmethod
    def update_potential_energy(obj):
        """
        Update/Load the potential_energy in the given obj from the attached storage

        Parameters
        ----------
        obj : Configuration
            the Configuration object to be updated

        """
        storage = obj._origin
        idx = obj.idx[storage]

        V = storage.variables['configuration_potential'][idx]
        potential_energy = u.Quantity(V, storage.units["configuration_potential"])

        obj.potential_energy = potential_energy

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

        return self.storage.variables['configuration_coordinates'][frame_indices,atom_indices,:].astype(np.float32).copy()

    def coordinates_as_array(self, frame_indices=None, atom_indices=None):
        '''
        Returns a numpy array consisting of all coordinates at the given indices

        Parameters
        ----------
        frame_indices : list of int
            configuration indices to be loaded
        atom_indices : list of int
            selects only the atoms to be returned. If None (Default) all atoms
            will be selected

        Returns
        -------
        numpy.ndarray, shape = (l,n)
            returns an array with `l` the number of frames and `n` the number
            of atoms
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
            selects only the atoms to be returned. If None (Default) all atoms
            will be selected


        Returns
        -------
        numpy.ndarray, shape = (l,n)
            returns an array with `l` the number of frames and `n` the number
            of atoms
        '''

        frame_indices = self.configuration_indices(idx)
        return self.coordinates_as_array(frame_indices, atom_indices)

    def _init(self):
        super(ConfigurationStore, self)._init()
        n_atoms = self.storage.n_atoms
        n_spatial = self.storage.n_spatial

        self.init_variable('configuration_coordinates', 'float',
                (self.db, 'atom','spatial'), self.dimension_units['length'],
                description="coordinates[configuration][atom][coordinate] " +
                            "are coordinate of atom 'atom' in dimension " +
                            "'coordinate' of configuration 'configuration'.",
                chunksizes=(1,n_atoms,n_spatial))

        self.init_variable('configuration_box_vectors', 'float',
                (self.db, 'spatial', 'spatial'),
                self.dimension_units['length'],
                chunksizes=(1,n_spatial,n_spatial))

        self.init_variable('configuration_potential', 'float',
                self.db,
                self.dimension_units['energy'],
                chunksizes=(1, ))