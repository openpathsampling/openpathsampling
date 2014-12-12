import numpy as np

from object_storage import ObjectStorage
from wrapper import savecache, loadcache
from opentis.trajectory import Trajectory


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
        return super(TrajectoryStorage, self).length('trajectory_snapshot', idx)

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
        begin = self.free_begin('trajectory_snapshot')
        length = len(trajectory)

        self.set_slice('trajectory', idx, begin, length)

        values = self.list_to_numpy(trajectory, 'snapshot')
        self.storage.variables['trajectory_snapshot_idx'][begin:begin+length] = values

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
        values = self.storage.variables['trajectory_snapshot_idx'][self.get_slice('trajectory', idx)]

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

        values = self.storage.variables['trajectory_snapshot_idx'][self.get_slice('trajectory', idx)]

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
        idx = storage.variables['trajectory_snapshot_idx'][:].astype(np.int32).copy()
        length = storage.variables['trajectory_snapshot_length'][:].astype(np.int32).copy()
        n_traj = self.count()

        return [ frames[idx[i]:idx[i] + length[i] ] for i in range(1, n_traj + 1) ]

    def _init(self, units=None):
        """
        Initialize the associated storage to allow for trajectory storage

        """
        super(TrajectoryStorage, self)._init()

        # index associated storage in class variable for all Trajectory instances to access
        ncfile = self.storage

        self.init_dimension('trajectory_snapshot')
        self.init_mixed('trajectory')
        self.init_variable('trajectory_snapshot_idx', 'index', 'trajectory_snapshot',
            description="trajectory[trajectory][frame] is the snapshot index (0..nspanshots-1) of frame 'frame' of trajectory 'trajectory'."
        )