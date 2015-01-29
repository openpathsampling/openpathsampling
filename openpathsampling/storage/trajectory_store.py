import numpy as np

from object_storage import ObjectStore
from openpathsampling.trajectory import Trajectory


class TrajectoryStore(ObjectStore):

    def __init__(self, storage):
        super(TrajectoryStore, self).__init__(storage, Trajectory)

    def save(self, trajectory, idx=None):
        """
        Add the current state of the trajectory to the storage. Saving also
        all referenced snapshots in it

        Parameters
        ----------
        trajectory : Trajectory()
            the trajectory to be saved
        idx : int or None
            if idx is not None the index will be used for saving in the storage.
            This might overwrite already existing trajectories!

        Notes
        -----
        This also saves all contained frames in the trajectory if not done yet.
        A single Trajectory object can only be saved once!
        """

        # Check if all snapshots are saved
        map(self.storage.snapshot.save, trajectory)

        values = self.list_to_numpy(trajectory, 'snapshot')
        self.storage.variables['trajectory_snapshot_idx'][idx] = values

    def snapshot_indices(self, idx):
        '''
        Load snapshot indices for trajectory with ID 'idx' from the storage

        Parameters
        ----------
        idx : int
            ID of the trajectory

        Returns
        -------
        list of int
            trajectory indices

        '''

        # get the values
        values = self.storage.variables['trajectory_snapshot_idx'][idx]

        # typecast to integer
        return self.list_from_numpy(values, 'index')

    def load(self, idx):
        '''
        Return a trajectory from the storage

        Parameters
        ----------
        idx : int
            index of the trajectory

        Returns
        -------
        Trajectory
            the trajectory

        '''

        values = self.storage.variables['trajectory_snapshot_idx'][idx]

        # typecast to snapshot
        snapshots = self.list_from_numpy(values, 'snapshot')

        trajectory = Trajectory(snapshots)

        return trajectory

    def iter_snapshot_indices(this, iter_range=None):
        '''
        Return an iterator over the lists of snapshot indices for all
        trajectories in the storage

        Parameters
        ----------
        iter_range : slice or None
            if this is not `None` it confines the iterator to objects specified
            in the slice

        Returns
        -------
        Iterator
            the iterator
        '''

        class ObjectIterator:
            def __init__(self):
                self.storage = this
                self.iter_range = iter_range
                if iter_range is None:
                    self.idx = 0
                    self.end = self.storage.count()
                else:
                    self.idx = iter_range.start
                    self.end = iter_range.stop

            def __iter__(self):
                return self

            def next(self):
                if self.idx < self.storage.count():
                    obj = self.snapshot_indices(self.idx)
                    if self.iter_range is not None and self.iter_range.step is not None:
                        self.idx += self.iter_range.step
                    else:
                        self.idx += 1
                    return obj
                else:
                    raise StopIteration()

        return ObjectIterator()

    def _init(self, units=None):
        super(TrajectoryStore, self)._init()

        # index associated storage in class variable for all Trajectory instances to access
        ncfile = self.storage

#        self.init_dimension('trajectory_snapshot')
#        self.init_mixed('trajectory')

        self.init_variable('trajectory_snapshot_idx', 'index', 'trajectory',
            description="trajectory[trajectory][frame] is the snapshot index (0..nspanshots-1) of frame 'frame' of trajectory 'trajectory'.",
            variable_length = True,
            chunksizes=(10240, )
        )