import abc
from uuid import UUID

from openpathsampling.netcdfplus import ObjectStore, LoaderProxy, StorableObject
import openpathsampling.engines as peng


# =============================================================================================
# ABSTRACT BASE CLASS FOR SNAPSHOTS
# =============================================================================================

class BaseSnapshotStore(ObjectStore):
    """
    An ObjectStore for Snapshots in netCDF files.
    """

    __metaclass__ = abc.ABCMeta

    def __init__(self, snapshot_class):
        """

        Attributes
        ----------
        snapshot_class : openpathsampling.BaseSnapshot
            a snapshot class that this Store is supposed to store

        """
        super(BaseSnapshotStore, self).__init__(peng.BaseSnapshot, json=False)
        self.snapshot_class = snapshot_class
        self._use_lazy_reversed = False
        if hasattr(snapshot_class, '__features__'):
            if '_reversed' in snapshot_class.__features__['lazy']:
                self._use_lazy_reversed = True

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
        snapshot : :obj:`BaseSnapshot`
            the loaded snapshot instance
        """

        # check if the reversed is in the cache
        try:
            return self.cache[BaseSnapshotStore.paired_idx(idx)].reversed
        except KeyError:
            pass

        # if not load and return it
        st_idx = int(idx / 2)

        obj = self.snapshot_class.__new__(self.snapshot_class)
        self.snapshot_class.init_empty(obj)

        self._get(st_idx, obj)
        if idx & 1:
            obj = obj.reversed

        # obj._reversed = LoaderProxy(self, BaseSnapshotStore.paired_idx(idx))
        return obj

    @abc.abstractmethod
    def _set(self, idx, snapshot):
        pass

    @abc.abstractmethod
    def _get(self, idx, snapshot):
        pass

    def _set_uuid(self, idx, uuid):
        if idx & 1:
            print 'Should not happen', idx
        self.storage.variables[self.prefix + '_uuid'][int(idx / 2)] = str(uuid)

    def _get_uuid(self, idx):
        uuid = UUID(self.storage.variables[self.prefix + '_uuid'][int(idx / 2)])
        if idx & 1:
            return UUID(int=int(uuid) ^ 1)
        else:
            return uuid

    def _update_uuid_in_cache(self, uuid, idx):
        # make sure to cast unicode to str
        uuid = str(uuid)
        if uuid != '':
            if uuid in self._uuid_idx:
                if self._uuid_idx[uuid] != int(idx):
                    raise RuntimeWarning('Already exists %s with idx %d -> %d' % (uuid, self._uuid_idx[uuid], int(idx)))
            else:
                self._uuid_idx[uuid] = int(idx)

            ruuid = str(UUID(int=int(UUID(uuid)) ^ 1))
            if ruuid in self._uuid_idx:
                if self._uuid_idx[ruuid] != int(idx) ^ 1:
                    raise RuntimeWarning('Inv Already exists %s with idx %d -> %d' % (ruuid, self._uuid_idx[ruuid], int(idx) ^ 1))
            else:
                self._uuid_idx[ruuid] = int(idx) ^ 1

    def update_uuid_cache(self):
        """
        Update the internal uuid cache with all stored uuids in the store.

        This allows to load by uuid for uuidd objects
        """
        if not self._uuids_loaded:
            for idx, uuid in enumerate(self.storage.variables[self.prefix + "_uuid"][:]):
                self._update_uuid_in_cache(uuid, idx * 2)

            self._uuids_loaded = True

    def _save(self, snapshot, idx):
        """
        Add the current state of the snapshot in the database.

        Parameters
        ----------
        snapshot :class:`openpathsampling.snapshots.AbstractSnapshot`
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
            # mark reversed as stored
            self.index[snapshot._reversed] = BaseSnapshotStore.paired_idx(idx)

    def save(self, obj, idx=None):
        if self.reference_by_uuid:
            ruuid = str(UUID(int=int(obj.__uuid__)))

            if ruuid in self.uuid_idx:
                # has been saved so quit and do nothing
                return obj.__uuid__

        return super(BaseSnapshotStore, self).save(obj, idx)

    def all(self):
        return peng.Trajectory(map(self.proxy, range(len(self))))

    def __len__(self):
        return 2 * super(BaseSnapshotStore, self).__len__()

    def duplicate(self, snapshot):
        """
        Store a duplicate of the snapshot as new

        Parameters
        ----------
        snapshot :class:`openpathsampling.snapshots.AbstractSnapshot`

        Returns
        -------
        int
            the index used for storing it in the store. This is the save as used by
            save.

        Notes
        -----
        This will circumvent the caching and indexing completely. This would be equivalent
        of creating a copy of the current snapshot and store this one and throw the copy
        away, leaving the given snapshot untouched. This allows you to treat the snapshot
        as mutual.

        The use becomes more obvious when applying to storing trajectories. The only way
        to make use of this feature is using the returned `idx`

        >>> idx = store.duplicate(snap)
        >>> loaded = store[idx]  # return a duplicated as new object
        >>> proxy = paths.LoaderProxy(store, idx) # use the duplicate without loading

        """
        idx = self.free()
        st_idx = int(idx / 2)
        self._set(st_idx, snapshot)

        return idx


# =============================================================================================
# FEATURE BASED SINGLE CLASS FOR ALL SNAPSHOT TYPES
# =============================================================================================

class FeatureSnapshotStore(BaseSnapshotStore):
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
