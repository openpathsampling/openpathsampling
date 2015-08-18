import numpy as np
from openpathsampling.snapshot import Snapshot, Configuration, Momentum
from object_storage import ObjectStore
from openpathsampling.trajectory import Trajectory
import simtk.unit as u

class SnapshotStore(ObjectStore):
    """
    An ObjectStore for Snapshots in netCDF files.
    """

    def __init__(self):
        super(SnapshotStore, self).__init__(Snapshot, json=False)

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

        configuration = self.vars['configuration'][idx]
        momentum = self.vars['momentum'][idx]
        momentum_reversed = self.vars['momentum_reversed'][idx]
        reversed_idx = self.vars['reversed_idx'][idx]

        snapshot = Snapshot(
            configuration=configuration,
            momentum=momentum,
            is_reversed=momentum_reversed,
            reversed_copy=None
        )
        snapshot_reversed = Snapshot(
            configuration=configuration,
            momentum=momentum,
            is_reversed=not momentum_reversed,
            reversed_copy=None
        )

        snapshot._reversed = snapshot_reversed
        snapshot_reversed._reversed = snapshot

        # fix caching!
        snapshot_reversed.idx[self.storage] = reversed_idx
        self.cache[reversed_idx] = snapshot_reversed

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

        self.vars['configuration'][idx] = snapshot.configuration
        self.vars['momentum'][idx] = snapshot.momentum
        self.vars['reversed'][idx] = snapshot._reversed

        self.vars['momentum_reversed'][idx] = snapshot.is_reversed

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
        return self.vars['configuration'][idx]

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
        return self.vars['momentum'][idx]

    def reversed_idx(self, idx):
        '''
        Load snapshot index for the reversed snapshot with ID 'idx'
        from the storage

        Parameters
        ----------
        idx : int
            index of the snapshot

        Returns
        -------
        int
            reversed snapshot indices
        '''
        return self.vars['reversed'][idx]

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
        return self.vars['momentum_reversed'][idx]


    def _init(self):
        '''
        Initializes the associated storage to index configuration_indices in it
        '''
        super(SnapshotStore, self)._init()

        self.init_variable('configuration', 'obj.configurations',
                description="the snapshot index (0..n_configuration-1) of snapshot '{idx}'.",
                chunksizes=(1, )
        )

        self.init_variable('momentum', 'obj.momenta',
                description="the snapshot index (0..n_momentum-1) 'frame' of snapshot '{idx}'.",
                chunksizes=(1, )
                )

        self.init_variable('momentum_reversed', 'bool', chunksizes=(1, ))

        self.init_variable('reversed', 'index',
                description="the idx of the reversed snapshot index (0..n_snapshot-1) 'snapshot' of snapshot '{idx}'.",
                chunksizes=(1, )
                )

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
            return obj.momenta.idx[self.storage]

        return idx


class MomentumStore(ObjectStore):
    """
    An ObjectStore for Momenta. Allows to store Momentum() instances in a netcdf file.
    """

    def __init__(self):
        super(MomentumStore, self).__init__(Momentum, json=False, load_partial=True)

        # attach delayed loaders
        self.set_variable_partial_loading('velocities')
        self.set_variable_partial_loading('kinetic_energy')

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
        if momentum.velocities is not None:
            self.vars['velocities'][idx,:,:] = momentum.velocities
        else:
            print 'ERROR : Momentum should not be empty'

        if momentum.kinetic_energy is not None:
            self.vars['kinetic_energy'][idx] = momentum.kinetic_energy
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


        velocities = self.vars['velocities'][idx]
        kinetic_energy = self.vars['kinetic_energy'][idx]

        momentum = Momentum(velocities=velocities, kinetic_energy=kinetic_energy)

        return momentum

    def load_empty(self, idx):
        momentum = Momentum()
        del momentum.velocities
        del momentum.kinetic_energy
        return momentum

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

        return self.variables['velocities'][frame_indices,atom_indices,:].astype(np.float32).copy()

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

        self.init_variable('velocities', 'float',
                dimensions=(self.prefix, 'atom', 'spatial'),
                units=self.dimension_units['velocity'],
                description="the velocity of atom 'atom' in dimension " +
                            "'coordinate' of momentum 'momentum'.",
                chunksizes=(1, n_atoms, n_spatial))

        self.init_variable('kinetic_energy', 'float', self.prefix,
                self.dimension_units['energy'],
                chunksizes=(1, ))


    
class ConfigurationStore(ObjectStore):
    def __init__(self):
        super(ConfigurationStore, self).__init__(Configuration, json=False, load_partial=True)

        # attach delayed loaders
        self.set_variable_partial_loading('coordinates')
        self.set_variable_partial_loading('box_vectors')
        self.set_variable_partial_loading('potential_energy')

    def save(self, configuration, idx = None):
        # Store configuration.
        self.vars['coordinates'][idx] = configuration.coordinates

        if configuration.potential_energy is not None:
            self.vars['potential'][idx] = configuration.potential_energy

        if configuration.box_vectors is not None:
            self.vars['box_vectors'][idx] = configuration.box_vectors

    def get(self, indices):
        return [ self.load(idx) for idx in indices ]

    def load(self, idx):
        coordinates = self.vars["coordinates"][idx]
        box_vectors = self.vars["box_vectors"][idx]
        potential_energy = self.vars["potential"][idx]

        configuration = Configuration(coordinates=coordinates, box_vectors = box_vectors, potential_energy=potential_energy)
        configuration.topology = self.storage.topology

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

        return self.storage.variables[self.prefix + '_coordinates'][frame_indices,atom_indices,:].astype(np.float32).copy()

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

        self.init_variable('coordinates', 'float',
                (self.prefix, 'atom','spatial'), self.dimension_units['length'],
                description="coordinate of atom '{ix[1]}' in dimension " +
                            "'{ix[2]}' of configuration '{ix[0]}'.",
                chunksizes=(1,n_atoms,n_spatial))

        self.init_variable('box_vectors', 'float',
                (self.prefix, 'spatial', 'spatial'),
                self.dimension_units['length'],
                chunksizes=(1,n_spatial,n_spatial))

        self.init_variable('potential', 'float',
                self.prefix,
                self.dimension_units['energy'],
                chunksizes=(1, ))