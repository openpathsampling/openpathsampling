from openpathsampling.engines.trajectory import Trajectory
from openpathsampling.netcdfplus import ObjectStore, LoaderProxy


class TrajectoryStore(ObjectStore):
    def __init__(self):
        super(TrajectoryStore, self).__init__(Trajectory)

    def to_dict(self):
        return {}

    def _save(self, trajectory, idx):
        self.vars['snapshots'][idx] = trajectory
        snapshot_store = self.storage.snapshots

        for frame, snapshot in enumerate(trajectory):
            if type(snapshot) is not LoaderProxy:
                loader = LoaderProxy(snapshot_store, snapshot_store.index[snapshot])
                trajectory[frame] = loader

    def _load(self, idx):
        trajectory = Trajectory(self.vars['snapshots'][idx])
        return trajectory

    def snapshot_indices(self, idx):
        """
        Load snapshot indices for trajectory with ID 'idx' from the storage

        Parameters
        ----------
        idx : int
            ID of the trajectory

        Returns
        -------
        list of int
            trajectory indices

        """

        # get the values
        return self.variables['snapshots'][idx].tolist()

    def iter_snapshot_indices(this, iter_range=None):
        """
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

        """

        class ObjectIterator:
            def __init__(self):
                self.storage = this
                self.iter_range = iter_range
                if iter_range is None:
                    self.idx = 0
                    self.end = len(self.storage)
                else:
                    self.idx = iter_range.start
                    self.end = iter_range.stop

            def __iter__(self):
                return self

            def next(self):
                if self.idx < len(self.storage):
                    obj = this.snapshot_indices(self.idx)
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

        self.create_variable('snapshots', 'lazyobj.snapshots',
                           dimensions=('...',),
                           description="trajectory[trajectory][frame] is the snapshot index " +
                                       "(0..nspanshots-1) of frame 'frame' of trajectory 'trajectory'.",
                           chunksizes=(10240,)
                           )
