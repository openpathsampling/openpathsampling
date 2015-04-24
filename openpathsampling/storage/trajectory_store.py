import numpy as np

from object_storage import ObjectStore
from openpathsampling.trajectory import Trajectory

# This adds delayed snapshot loading support to Trajectory

def load_missing_snapshot(func):
    def getter(self, *args, **kwargs):
        item = func(self, *args, **kwargs)

        if type(item) is tuple:
            item = item[0][item[1]]

        return item

    return getter

Trajectory.__getitem__ = load_missing_snapshot(Trajectory.__getitem__)
Trajectory.__getslice__ = load_missing_snapshot(Trajectory.__getslice__)

class TrajectoryStore(ObjectStore):
    def __init__(self, storage, lazy=True, use_snapshot_cache=True):
        super(TrajectoryStore, self).__init__(storage, Trajectory)
        self.lazy = lazy
        self.use_snapshot_cache = use_snapshot_cache

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
        values = []
        for snap in list.__iter__(trajectory):
            if type(snap) is not tuple:
                self.storage.snapshots.save(snap)
                values.append(snap.idx[self.storage])
            elif snap[0].storage is not self.storage:
                new_snap = snap[0][snap[1]]
                self.storage.snapshots.save(new_snap)
                values.append(new_snap.idx[self.storage])
            else:
                values.append(snap[1])

#        map(self.storage.snapshots.save, trajectory)
#        values = self.list_to_numpy(trajectory, 'snapshot')

        values = self.list_to_numpy(values, 'index')
        self.storage.variables['trajectory_snapshot_idx'][idx] = values

        # self.storage.sync()

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
        if self.lazy:
            if self.use_snapshot_cache:
                snapshots = [ self.storage.snapshots.cache[idx] if idx in self.storage.snapshots.cache else tuple([self.storage.snapshots, idx]) for idx in self.list_from_numpy(values, 'int') ]
            else:
                snapshots = [ tuple([self.storage.snapshots, idx]) for idx in self.list_from_numpy(values, 'int') ]
        else:
            snapshots = self.list_from_numpy(values, 'snapshots')

        trajectory = Trajectory(snapshots)
        # save the used storage to load snapshots if required

        trajectory.storage = self.storage

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

        self.init_variable('trajectory_snapshot_idx', 'index', 'trajectory',
            description="trajectory[trajectory][frame] is the snapshot index (0..nspanshots-1) of frame 'frame' of trajectory 'trajectory'.",
            variable_length = True,
            chunksizes=(10240, )
        )