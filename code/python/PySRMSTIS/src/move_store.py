import numpy as np
import mdtraj as md

from object_storage import ObjectStorage
from move import Move

class MoveStorage(ObjectStorage):

    def __init__(self, storage):
        super(MoveStorage, self).__init__(storage, Move)

    def length(self, idx):
        '''
        Return the length of a move from the storage

        Parameters
        ----------
        idx : int
            index of the move

        Returns
        -------
        length : int
            Number of frames in the move
        '''
        return int(self.storage.variables['move_length'][idx])

    def save(self, move, idx=None):
        """
        Add the current state of the move in the database. If nothing has changed then the move gets stored using the same snapshots as before. Saving lots of diskspace

        Parameters
        ----------
        move : Move()
            the move to be saved
        idx : int or None
            if idx is not None the index will be used for saving in the storage. This might overwrite already existing trajectories!

        Notes
        -----
        This also saves all contained frames in the move if not done yet.
        A single Move object can only be saved once!
        """

        idx = super(MoveStorage, self).index(move, idx)

        if idx is not None:
            storage = self.storage

            begin = self._free_idx()

    #        print 'Begin :', begin
    #        print 'Index :', idx

            nframes = len(move)
            for frame_index in range(nframes):
                frame = move[frame_index]
                storage.snapshot.save(frame)

#                print 'Position :', begin + frame_index

                storage.variables['move_configuration_idx'][begin + frame_index] = frame.configuration.idx[storage]
                storage.variables['move_momentum_idx'][begin + frame_index] = frame.momentum.idx[storage]
                storage.variables['move_momentum_reversed'][begin + frame_index] = frame.reversed

            storage.variables['move_length'][idx] = nframes
            storage.variables['move_idx'][idx] = begin

        return

    def load(self, idx, momentum = True):
        '''
        Return a move from the storage

        Parameters
        ----------
        idx : int
            index of the move (counts from 1)

        Returns
        -------
        move : Move
            the move
        '''
        frames_c =  self.configuration_indices(idx)
        frames_m =  self.momentum_indices(idx)
        reversed_m =  self.momentum_reversed(idx)
        obj = self.from_indices(frames_c, frames_m, reversed_m)
        obj.idx[self.storage] = idx

        return obj



    def from_indices(self, frames_configuration, frames_momentum, momenta_reversed, storage = None):
        '''
        Return a move from the storage constructed from a list of snapshot indices

        Parameters
        ----------
        frames_configuration : list of int
            list of configuration indices to be used to generate the move
        frames_momentum : list of int
            list of momentum indices to be used to generate the move
        momentum_reversed : list of bool
            list indicating if frames are reversed

        Returns
        -------
        move : Move
            the move
        '''
        move = Move()

        if frames_momentum is not None:
            for idcs in zip(frames_configuration, frames_momentum, momenta_reversed):
                snapshot = self.storage.snapshot.load(*idcs)
                move.append(snapshot)
        else:
            for frame in zip(frames_configuration, momenta_reversed):
                snapshot = storage.snapshot( frame[0], None , frame[2])
                move.append(snapshot)

        return move

    def _free_idx(self):
        '''
        Return the number of the next free ID

        Returns
        -------
        index : int
            the number of the next free index in the storage. Used to store a new snapshot.
        '''
        length = int(len(self.storage.dimensions['frames']))
        return length + 1

    def all_momentum_indices(self):
        '''
        Return a list of frame indices for all trajectories in the storage

        Returns
        -------
        list : list of list of int
            a list of list of frame IDs for all trajectories.
        '''


        storage = self.storage
        frames = storage.variables['move_momentum_idx'][:].astype(np.int32).copy()
        idx = storage.variables['move_length'][:].astype(np.int32).copy()
        length = storage.variables['move_length'][:].astype(np.int32).copy()
        n_traj =  self.count(storage)

        return [ frames[idx[i]:idx[i] + length[i] ] for i in range(1, n_traj + 1) ]


    def all_configuration_indices(self):
        '''
        Return a list of frame indices for all trajectories in the storage

        Returns
        -------
        list : list of list of int
            a list of list of frame IDs for all trajectories.
        '''

        storage = self.storage
        frames = storage.variables['move_configuration_idx'][:].astype(np.int32).copy()
        idx = storage.variables['move_length'][:].astype(np.int32).copy()
        length = storage.variables['move_length'][:].astype(np.int32).copy()
        n_traj =  self.count(storage)

        return [ frames[idx[i]:idx[i] + length[i] ] for i in range(1, n_traj + 1) ]

    def all_snapshot_coordinates_as_mdtraj(self, atom_indices = None):
        """
        Return all snapshots as a mdtraj.Move object using only the specified atoms

        Parameters
        ----------
        atom_indices (list of int, Default: None) - list of atom indices to be used for the move

        """

        output = self.coordinates_as_array(atom_indices = atom_indices)

        topology = self.storage.topology

        if atom_indices is not None:
            topology = topology.subset(atom_indices)

        return md.Move(output, topology)

    def velocities_as_array(self, idx, atom_indices=None):
        '''
        Returns a numpy array consisting of all velocities of a move

        Parameters
        ----------
        idx : int
            index of the move to be loaded
        atom_indices : list of int
            selects only the atoms to be returned. If None (Default) all atoms will be selected

        Returns
        -------
        numpy.ndarray, shape = (l,n)
            returns an array with `l` the number of frames and `n` the number of atoms
        '''

        frame_indices = self.configuration_indices(idx)
        return self.storage.momentum.velocities_as_array(frame_indices, atom_indices)


    def _init(self):
        """
        Initialize the associated storage to allow for move storage

        """
        super(MoveStorage, self)._init()

        # index associated storage in class variable for all Move instances to access
        ncfile = self.storage

        # define dimensions used in trajectories
        ncfile.createDimension('frames', 0)                     # unlimited number of iterations

        # Create variables for trajectories
        ncvar_move_configuration_idx  = ncfile.createVariable('move_configuration_idx', 'u4', 'frames')
        ncvar_move_momentum_idx       = ncfile.createVariable('move_momentum_idx', 'u4', 'frames')
        ncvar_move_momentum_reversed  = ncfile.createVariable('move_momentum_reversed', 'b', 'frames')
        ncvar_move_path_hamiltonian   = ncfile.createVariable('move_path_hamiltonians', 'f', 'frames')
        ncvar_move_length             = ncfile.createVariable('move_length', 'u4', self.idx_dimension)
        ncvar_move_idx                = ncfile.createVariable('move_idx', 'u4', self.idx_dimension)


        # Define units for snapshot variables.
        setattr(ncvar_move_path_hamiltonian,      'units', 'none')
        setattr(ncvar_move_configuration_idx,     'units', 'none')
        setattr(ncvar_move_momentum_idx,          'units', 'none')
        setattr(ncvar_move_length,                'units', 'none')
        setattr(ncvar_move_idx,                   'units', 'none')
        setattr(ncvar_move_momentum_reversed,     'units', 'none')

        # Define long (human-readable) names for variables.
        setattr(ncvar_move_configuration_idx,    "long_name", "move[move][frame] is the snapshot index (0..nspanshots-1) of frame 'frame' of move 'move'.")
        setattr(ncvar_move_momentum_idx,         "long_name", "move[move][frame] is the snapshot index (0..nspanshots-1) of frame 'frame' of move 'move'.")
