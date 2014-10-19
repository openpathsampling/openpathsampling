import numpy as np
import mdtraj as md

from object_storage import ObjectStorage
from trajectory import Trajectory

class TrajectoryStorage(ObjectStorage):

    def __init__(self, storage):
        super(TrajectoryStorage, self).__init__(storage, Trajectory)

    def _slice(self, idx):
        idx = self._idx(idx)
        end =  idx + self.length(idx)
        return slice(idx,end)

    def _idx(self, idx):
        return int(self.storage.variables['trajectory_idx'][idx])

    def length(self, idx):
        '''
        Return the length of a trajectory from the storage

        Parameters
        ----------
        idx : int
            index of the trajectory

        Returns
        -------
        length : int
            Number of frames in the trajectory
        '''
        return int(self.storage.variables['trajectory_length'][idx])

    def index(self, trajectory, idx=None):
        """
        Add the current state of the trajectory in the database. If nothing has changed then the trajectory gets stored using the same snapshots as before. Saving lots of diskspace

        Parameters
        ----------
        trajectory : Trajectory()
            the trajectory to be saved
        idx : int or None
            if idx is not None the index will be used for saving in the storage. This might overwrite already existing trajectories!

        Notes
        -----
        This also saves all contained frames in the trajectory if not done yet.
        A single Trajectory object can only be saved once!
        """

        idx = super(TrajectoryStorage, self).index(trajectory, idx)

        if idx is not None:
            storage = self.storage

            begin = self._free_idx()

    #        print 'Begin :', _idx
    #        print 'Index :', idx

            nframes = len(trajectory)
            for frame_index in range(nframes):
                frame = trajectory[frame_index]
                storage.snapshot.save(frame)

#                print 'Position :', begin + frame_index

                storage.variables['trajectory_configuration_idx'][begin + frame_index] = frame.configuration.idx[storage]
                storage.variables['trajectory_momentum_idx'][begin + frame_index] = frame.momentum.idx[storage]
                storage.variables['trajectory_momentum_reversed'][begin + frame_index] = frame.reversed

            storage.variables['trajectory_length'][idx] = nframes
            storage.variables['trajectory_idx'][idx] = begin

        return


    def momentum_indices(self, idx):
        '''
        Load trajectory indices for trajectory with ID 'idx' from the storage

        ARGUMENTS

        idx (int) - ID of the trajectory

        Returns
        -------
        trajectory (list of int) - trajectory indices
        '''
        return self.storage.variables['trajectory_momentum_idx'][self._slice(idx)]



    def momentum_reversed(self, idx):
        '''
        Load trajectory with ID 'idx' from the storage and return a list of reversed indicators for the momenta

        Parameters
        ----------

        idx : int
            index of the trajectory

        Returns
        -------
        list of boolean
            list of boolean which frames in the trajectory are reversed
        '''
        return self.storage.variables['trajectory_momentum_reversed'][self._slice(idx)]


    def configuration_indices(self, idx):
        '''
        Load trajectory indices for trajectory with ID 'idx' from the storage

        ARGUMENTS

        idx (int) - ID of the trajectory

        Returns
        -------
        trajectory (list of int) - trajectory indices
        '''
        return self.storage.variables['trajectory_configuration_idx'][self._slice(idx)]



    def indices(self, idx):
        '''
        Load trajectory indices for trajectory with ID 'idx' from the storage

        ARGUMENTS

        idx (int) - ID of the trajectory

        Returns
        -------
        trajectory (list of int) - trajectory indices
        '''

        return (
                self.configuration_indices(idx),
                self.momentum_indices(idx)
        )


    def load(self, idx, momentum = True):
        '''
        Return a trajectory from the storage

        Parameters
        ----------
        idx : int
            index of the trajectory (counts from 1)

        Returns
        -------
        trajectory : Trajectory
            the trajectory
        '''
        frames_c =  self.configuration_indices(idx)
        frames_m =  self.momentum_indices(idx)
        reversed_m =  self.momentum_reversed(idx)
        obj = self.from_indices(frames_c, frames_m, reversed_m)
        obj.idx[self.storage] = idx

        return obj



    def from_indices(self, frames_configuration, frames_momentum, momenta_reversed, storage = None):
        '''
        Return a trajectory from the storage constructed from a list of snapshot indices

        Parameters
        ----------
        frames_configuration : list of int
            list of configuration indices to be used to generate the trajectory
        frames_momentum : list of int
            list of momentum indices to be used to generate the trajectory
        momentum_reversed : list of bool
            list indicating if frames are reversed

        Returns
        -------
        trajectory : Trajectory
            the trajectory
        '''
        trajectory = Trajectory()

        if frames_momentum is not None:
            for idcs in zip(frames_configuration, frames_momentum, momenta_reversed):
                snapshot = self.storage.snapshot.load(*idcs)
                trajectory.append(snapshot)
        else:
            for frame in zip(frames_configuration, momenta_reversed):
                snapshot = storage.snapshot( frame[0], None , frame[2])
                trajectory.append(snapshot)

        return trajectory

    def _free_idx(self):
        '''
        Return the number of the next _free_idx ID

        Returns
        -------
        index : int
            the number of the next _free_idx index in the storage. Used to store a new snapshot.
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
        frames = storage.variables['trajectory_momentum_idx'][:].astype(np.int32).copy()
        idx = storage.variables['trajectory_length'][:].astype(np.int32).copy()
        length = storage.variables['trajectory_length'][:].astype(np.int32).copy()
        n_traj =  self.number(storage)

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
        frames = storage.variables['trajectory_configuration_idx'][:].astype(np.int32).copy()
        idx = storage.variables['trajectory_length'][:].astype(np.int32).copy()
        length = storage.variables['trajectory_length'][:].astype(np.int32).copy()
        n_traj =  self.number(storage)

        return [ frames[idx[i]:idx[i] + length[i] ] for i in range(1, n_traj + 1) ]

    def all_snapshot_coordinates_as_mdtraj(self, atom_indices = None):
        """
        Return all snapshots as a mdtraj.Trajectory object using only the specified atoms

        Parameters
        ----------
        atom_indices (list of int, Default: None) - list of atom indices to be used for the trajectory

        """

        output = self.coordinates_as_array(atom_indices = atom_indices)

        topology = self.storage.topology

        if atom_indices is not None:
            topology = topology.subset(atom_indices)

        return md.Trajectory(output, topology)

    def velocities_as_array(self, idx, atom_indices=None):
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

        frame_indices = self.configuration_indices(idx)
        return self.storage.momentum.velocities_as_array(frame_indices, atom_indices)


    def _init(self):
        """
        Initialize the associated storage to allow for trajectory storage

        """
        super(TrajectoryStorage, self)._init()

        # index associated storage in class variable for all Trajectory instances to access
        ncfile = self.storage

        # define dimensions used in trajectories
        ncfile.createDimension('frames', 0)                     # unlimited number of iterations

        # Create variables for trajectories
        ncvar_trajectory_configuration_idx  = ncfile.createVariable('trajectory_configuration_idx', 'u4', 'frames')
        ncvar_trajectory_momentum_idx       = ncfile.createVariable('trajectory_momentum_idx', 'u4', 'frames')
        ncvar_trajectory_momentum_reversed  = ncfile.createVariable('trajectory_momentum_reversed', 'b', 'frames')
        ncvar_trajectory_path_hamiltonian   = ncfile.createVariable('trajectory_path_hamiltonians', 'f', 'frames')
        ncvar_trajectory_length             = ncfile.createVariable('trajectory_length', 'u4', self.idx_dimension)
        ncvar_trajectory_idx                = ncfile.createVariable('trajectory_idx', 'u4', self.idx_dimension)


        # Define units for snapshot variables.
        setattr(ncvar_trajectory_path_hamiltonian,      'units', 'none')
        setattr(ncvar_trajectory_configuration_idx,     'units', 'none')
        setattr(ncvar_trajectory_momentum_idx,          'units', 'none')
        setattr(ncvar_trajectory_length,                'units', 'none')
        setattr(ncvar_trajectory_idx,                   'units', 'none')
        setattr(ncvar_trajectory_momentum_reversed,     'units', 'none')

        # Define long (human-readable) names for variables.
        setattr(ncvar_trajectory_configuration_idx,    "long_name", "trajectory[trajectory][frame] is the snapshot index (0..nspanshots-1) of frame 'frame' of trajectory 'trajectory'.")
        setattr(ncvar_trajectory_momentum_idx,         "long_name", "trajectory[trajectory][frame] is the snapshot index (0..nspanshots-1) of frame 'frame' of trajectory 'trajectory'.")
