import numpy as np
import mdtraj as md

from object_storage import ObjectStorage
from wrapper import savecache, loadcache
from trajectory import Trajectory, Sample
from snapshot import Configuration, Momentum, Snapshot


class TrajectoryStorage(ObjectStorage):

    def __init__(self, storage):
        super(TrajectoryStorage, self).__init__(storage, Trajectory)

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
        return int(self.storage.variables['trajectory_snapshots_length'][idx])

    @savecache
    def save(self, trajectory, idx=None):
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

        # Check if all snapshots are saved
        map(self.storage.snapshot.save, trajectory)

        # Find a free position to save snapshot ids
        begin = self.free_begin('trajectory_snapshots')
        length = len(trajectory)

        self.set_slice('trajectory_snapshots', idx, begin, length)
        self.set_list_as_type('trajectory_snapshots', idx, begin, length, 'snapshot')

    def snapshot_indices(self, idx):
        '''
        Load snapshot indices for trajectory with ID 'idx' from the storage

        ARGUMENTS

        idx (int) - ID of the trajectory

        Returns
        -------
        list of int - trajectory indices
        '''

        # get the values
        values = self.storage.variables['trajectory_snapshots'][idx, self.get_slice('trajectory_snapshots', idx)]

        # typecast to integer
        return self.list_from_numpy(values, 'index')

    @loadcache
    def load(self, idx, lazy = None):
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

        values = self.storage.variables['trajectory_snapshots'][idx, self.get_slice('trajectory_snapshots', idx)]

        # typecast to snapshot
        snapshots = self.list_from_numpy(values, 'snapshot')

        trajectory = Trajectory(snapshots)
        trajectory.idx[self.storage] = idx

        return trajectory

    def all_snapshot_indices(self):
        '''
        Return a list of snapshot indices for all trajectories in the storage

        Returns
        -------
        list : list of list of int
            a list of list of frame IDs for all trajectories.
        '''

        storage = self.storage
        frames = storage.variables['trajectory_snapshot_idx'][:].astype(np.int32).copy().tolist()
        idx = storage.variables['trajectory_snapshots_idx'][:].astype(np.int32).copy()
        length = storage.variables['trajectory_snapshots_length'][:].astype(np.int32).copy()
        n_traj = self.count()

        return [ frames[idx[i]:idx[i] + length[i] ] for i in range(1, n_traj + 1) ]

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

        frame_indices = self.momentum_indices(idx)
        return self.storage.momentum.velocities_as_array(frame_indices, atom_indices)


    def _init(self):
        """
        Initialize the associated storage to allow for trajectory storage

        """
        super(TrajectoryStorage, self)._init()

        # index associated storage in class variable for all Trajectory instances to access
        ncfile = self.storage

        self.init_dimension('trajectory_snapshots')
        self.init_variable('trajectory_snapshot', 'index', 'trajectory_snapshots',
            description="trajectory[trajectory][frame] is the snapshot index (0..nspanshots-1) of frame 'frame' of trajectory 'trajectory'."
        )

class SampleStorage(ObjectStorage):
    def __init__(self, storage):
        super(SampleStorage, self).__init__(storage, Sample)

    @savecache
    def save(self, origin, idx=None):
        """
        Add the current state of the origin in the database. If nothing has changed then the origin gets stored using the same snapshots as before. Saving lots of diskspace

        Parameters
        ----------
        origin : Sample()
            the origin to be saved
        idx : int or None
            if idx is not None the index will be used for saving in the storage. This might overwrite already existing trajectories!

        Notes
        -----
        This also saves all contained frames in the origin if not done yet.
        A single Sample object can only be saved once!
        """

        if idx is not None:
            self.storage.trajectory.save(origin.trajectory)
            self.set_object('origin_trajectory', idx, origin.trajectory)

            self.storage.ensemble.save(origin.ensemble)
            self.set_object('origin_ensemble', idx, origin.ensemble)

            self.storage.pathmover.save(origin.mover)
            self.set_object('origin_mover', idx, origin.mover)

            self.storage.movedetails.save(origin.details)
            self.set_object('origin_details', idx, origin.details)

    @loadcache
    def load(self, idx, momentum = True):
        '''
        Return a origin from the storage

        Parameters
        ----------
        idx : int
            index of the origin (counts from 1)

        Returns
        -------
        origin : Sample
            the origin
        '''
        trajectory_idx = int(self.storage.variables['origin_trajectory_idx'][idx])
        ensemble_idx = int(self.storage.variables['origin_ensemble_idx'][idx])
        mover_idx = int(self.storage.variables['origin_mover_idx'][idx])
        details_idx = int(self.storage.variables['origin_details_idx'][idx])

        obj = Sample(
            trajectory=self.storage.trajectory.load(trajectory_idx, lazy=True),
            mover=self.storage.pathmover.load(mover_idx, lazy=True),
            ensemble=self.storage.ensemble.load(ensemble_idx),
            details=self.storage.movedetails.load(details_idx)
        )

        return obj

    def _init(self):
        """
        Initialize the associated storage to allow for origin storage

        """
        super(SampleStorage, self)._init()

        # New short-hand definition
        self.init_variable('origin_trajectory_idx', 'u4')
        self.init_variable('origin_ensemble_idx', 'u4')
        self.init_variable('origin_mover_idx', 'u4')
        self.init_variable('origin_details_idx', 'u4')