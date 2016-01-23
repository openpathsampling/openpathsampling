from openpathsampling.netcdfplus import ObjectStore, LoaderProxy
from openpathsampling.snapshot import Snapshot, AbstractSnapshot
from openpathsampling.trajectory import Trajectory

import abc


# =============================================================================================
# ABSTRACT BASE CLASS FOR SNAPSHOTS
# =============================================================================================

class AbstractSnapshotStore(ObjectStore):
    """
    An ObjectStore for Snapshots in netCDF files.
    """

    __metaclass__ = abc.ABCMeta

    def __init__(self, snapshot_class):
        """

        Attributes
        ----------
        snapshot_class : openpathsampling.AbstractSnapshot
            a snapshot class that this Store is supposed to store

        """
        super(AbstractSnapshotStore, self).__init__(AbstractSnapshot, json=False)
        self.snapshot_class = snapshot_class

    def __repr__(self):
        return "store.%s[%s(%s)]" % (
            self.prefix, self.snapshot_class.__name__, self.content_class.__name__)

    @staticmethod
    def paired_idx(idx):
        """
        Return the paired index

        Snapshots are stored in pairs (2n, 2n+1) where one is the reversed copy.
        This make storing CVs easier. This function allows to get the paired index
        or the index of snapshot.reversed

        The implementation uses the trick that all you have to do is flip the lowest bit
        that determines even or odd.

        Parameters
        ----------
        idx : int
            the one part of the paired index

        Returns
        -------
        int
            the other part of the paired index
        """
        return idx ^ 1

    def to_dict(self):
        return {
            'snapshot_class': self.snapshot_class
        }

    def _load(self, idx):
        """
        Load a snapshot from the storage.

        Parameters
        ----------
        idx : int
            the integer index of the snapshot to be loaded

        Returns
        -------
        snapshot : Snapshot
            the loaded snapshot instance
        """

        # check if the reversed is in the cache
        try:
            obj = self.cache[AbstractSnapshotStore.paired_idx(idx)].create_reversed()
            obj._reversed = LoaderProxy(self, AbstractSnapshotStore.paired_idx(idx))
            return obj
        except KeyError:
            pass

        # if not load and return it
        st_idx = int(idx / 2)

        obj = self.snapshot_class.__new__(self.snapshot_class)
        AbstractSnapshot.__init__(obj)

        self._get(st_idx, obj)
        if idx & 1:
            obj = obj.create_reversed()

        obj._reversed = LoaderProxy(self, AbstractSnapshotStore.paired_idx(idx))
        return obj


    @abc.abstractmethod
    def _set(self, idx, snapshot):
        pass

    @abc.abstractmethod
    def _get(self, idx, snapshot):
        pass

    def _save(self, snapshot, idx):
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


        st_idx = int(idx / 2)

        self._set(st_idx, snapshot)

        if snapshot._reversed is not None:
            snapshot._reversed._reversed = LoaderProxy(self, idx)
            # mark reversed as stored
            self.index[snapshot._reversed] = AbstractSnapshotStore.paired_idx(idx)

        snapshot._reversed = LoaderProxy(self, AbstractSnapshotStore.paired_idx(idx))

    def all(self):
        return Trajectory([LoaderProxy(self, idx) for idx in range(len(self))])

    def __len__(self):
        return 2 * super(AbstractSnapshotStore, self).__len__()

    def duplicate(self, snapshot):
        """
        Store a duplicate of the snapshot as new

        Parameters
        ----------
        snapshot

        Returns
        -------

        """
        idx = self.free()
        st_idx = int(idx / 2)
        self._set(st_idx, snapshot.create_reversed())

        return idx



# =============================================================================================
# FEATURE BASED SINGLE CLASS FOR ALL SNAPSHOT TYPES
# =============================================================================================

class FeatureSnapshotStore(AbstractSnapshotStore):
    """
    An ObjectStore for Snapshots in netCDF files.
    """

    def __init__(self, snapshot_class):
        super(FeatureSnapshotStore, self).__init__(snapshot_class)

    @property
    def classes(self):
        return self.snapshot_class.__features__['classes']

    @property
    def parameters(self):
        return self.snapshot_class.__features__['parameters']

    def _set(self, idx, snapshot):
        [self.write(attr, idx, snapshot) for attr in self.parameters]

    def _get(self, idx, snapshot):
        [setattr(snapshot, attr, self.vars[attr][idx]) for attr in self.parameters]

    def _init(self):
        super(FeatureSnapshotStore, self)._init()

        for feature in self.classes:
            if hasattr(feature, 'netcdfplus_init'):
                feature.netcdfplus_init(self)
