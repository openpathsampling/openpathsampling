from openpathsampling.engines.trajectory import Trajectory
from openpathsampling.netcdfplus import ObjectStore, LoaderProxy


class TrajectoryStore(ObjectStore):
    def __init__(self):
        super(TrajectoryStore, self).__init__(Trajectory)

    def to_dict(self):
        return {}

    def _save(self, trajectory, idx):
        self.vars['snapshots'][idx] = trajectory
        store = self.storage.snapshots

        for frame, snapshot in enumerate(trajectory.iter_proxies()):
            if type(snapshot) is not LoaderProxy:
                loader = store.proxy(snapshot)
                trajectory[frame] = loader

    def mention(self, trajectory):
        """
        Save a trajectory and store its snapshots only shallow

        This will mention the ids of all snapshots in the file but not save
        the content of all the snapshots. This way you can store CV values
        if you want

        Parameters
        ----------
        trajectory : :class:`openpathsampling.Trajectory`

        """
        snap_store = self.storage.snapshots
        current_mention = snap_store.only_mention
        snap_store.only_mention = True
        self.save(trajectory)
        snap_store.only_mention = current_mention

    def _load(self, idx):
        trajectory = Trajectory(self.vars['snapshots'][idx])
        return trajectory

    def cache_all(self):
        """Load all samples as fast as possible into the cache

        """
        if not self._cached_all:
            idxs = range(len(self))
            snaps = self.vars['snapshots'][:]

            [self.add_single_to_cache(i, j) for i, j in zip(
                idxs,
                snaps)]

            self._cached_all = True

    def add_single_to_cache(self, idx, snaps):
        """
        Add a single object to cache by json

        Parameters
        ----------
        idx : int
            the index where the object was stored
        snaps : list of `BaseSnapshot`
            json string the represents a serialized version of the stored object
        """

        if idx not in self.cache:
            obj = Trajectory(snaps)

            self._get_id(idx, obj)

            self.cache[idx] = obj
            self.index[obj.__uuid__] = idx

            return obj

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

    def iter_snapshot_indices(self):
        """
        Return an iterator over the lists of snapshot indices for all
        trajectories in the storage

        Returns
        -------
        Iterator
            the iterator

        """
        for snap_idx in range(len(self)):
            yield self.snapshot_indices(snap_idx)

    def initialize(self, units=None):
        super(TrajectoryStore, self).initialize()

        # index associated storage in class variable for all Trajectory
        # instances to access

        self.create_variable(
            'snapshots',
            'lazyobj.snapshots',
            dimensions=('...',),
            description="trajectory[trajectory][frame] is the snapshot index "
                        "(0..nspanshots-1) of frame 'frame' of trajectory "
                        "'trajectory'.",
            chunksizes=(65536,)
        )
