import abc
from uuid import UUID

from openpathsampling.netcdfplus import ObjectStore, LoaderProxy, StorableObject
from openpathsampling.netcdfplus.objects import UUIDDict, IndexedObjectStore
import openpathsampling.engines as peng

from collections import OrderedDict


# =============================================================================================
# ABSTRACT BASE CLASS FOR SNAPSHOTS
# =============================================================================================

class UUIDReversalDict(UUIDDict):
    @staticmethod
    def rev_id(obj):
        return StorableObject.ruuid(UUIDReversalDict.id(obj))

    def __setitem__(self, key, value):
        OrderedDict.__setitem__(self, self.id(key), value)
        OrderedDict.__setitem__(self, self.rev_id(key), value ^ 1)

    def __delitem__(self, key):
        OrderedDict.__delitem__(self, self.id(key))
        OrderedDict.__delitem__(self, self.rev_id(key))


class BaseSnapshotStore(IndexedObjectStore):
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
            if '_reversed' in snapshot_class.__features__.lazy:
                self._use_lazy_reversed = True

    def create_uuid_index(self):
        return UUIDReversalDict()

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

    def _set_id(self, idx, obj):
        if self.reference_by_uuid:
            self.vars['uuid'][int(idx / 2)] = obj.__uuid__

    def _get_id(self, idx, obj):
        if self.reference_by_uuid:
            uuid = self.vars['uuid'][int(idx / 2)]
            if idx & 1:
                uuid = StorableObject.ruuid(uuid)

            obj.__uuid__ = uuid

    def load_indices(self):
        if self.reference_by_uuid:
            for idx, uuid in enumerate(self.vars['uuid'][:]):
                self.index[uuid] = idx * 2

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

        if snapshot._reversed is not None:
            if not self.reference_by_uuid and snapshot._reversed in self.index:
                # seems we have already stored this snapshot but didn't know about it
                raise RuntimeWarning('This should never happen! Please report a bug!')
            else:
                # mark reversed as stored
                self.index[snapshot._reversed] = BaseSnapshotStore.paired_idx(idx)

        self._set(st_idx, snapshot)

        if snapshot._reversed is not None:
            # mark reversed as stored
            self.index[snapshot._reversed] = BaseSnapshotStore.paired_idx(idx)

    def save(self, obj, idx=None):
        if self.reference_by_uuid:
            ruuid = str(UUID(int=int(obj.__uuid__)))

            if ruuid in self.index:
                # has been saved so quit and do nothing
                return obj.__uuid__

        if obj._reversed is not None:
            if not self.reference_by_uuid and obj._reversed in self.index:
                # the reversed copy has been saved so quit and return the paired idx
                self.index[obj] = BaseSnapshotStore.paired_idx(self.index[obj._reversed])

        return super(BaseSnapshotStore, self).save(obj, idx)

    def all(self):
        if self.reference_by_uuid:
            return peng.Trajectory(map(self.proxy, self.index))
        else:
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

    def idx(self, obj):
        try:
            return self.index[obj]
        except KeyError:
            pass

        try:
            return BaseSnapshotStore.paired_idx(self.index[obj.reversed])
        except KeyError:
            return None


class BaseSnapshotIndexedStore(IndexedObjectStore):
    """
    An ObjectStore for Snapshots in netCDF files.
    """

    __metaclass__ = abc.ABCMeta

    def __init__(self, descriptor):
        """

        Attributes
        ----------
        snapshot_class : openpathsampling.BaseSnapshot
            a snapshot class that this Store is supposed to store

        """
        super(BaseSnapshotIndexedStore, self).__init__(peng.BaseSnapshot, json=False)
        self.descriptor = descriptor
        self._dimensions = descriptor.dimensions
        self._cls = self.descriptor.snapshot_class

        # self._use_lazy_reversed = False
        # if hasattr(snapshot_class, '__features__'):
        #     if '_reversed' in snapshot_class.__features__.lazy:
        #         self._use_lazy_reversed = True

    @property
    def reference_by_uuid(self):
        # This one does explicitly use integer indices
        return False

    @property
    def snapshot_class(self):
        return self._cls

    @property
    def dimensions(self):
        return self.descriptor['dimensions']

    def __repr__(self):
        return "store.%s[%s(%s)]" % (
            self.prefix, self._cls.__name__, self.content_class.__name__)

    def to_dict(self):
        return {
            'descriptor': self.descriptor,
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

        # if not load and return it
        st_idx = int(idx)

        obj = self._cls.__new__(self._cls)
        self._cls.init_empty(obj)
        self._get(st_idx, obj)
        return obj

    @abc.abstractmethod
    def _set(self, idx, snapshot):
        pass

    @abc.abstractmethod
    def _get(self, idx, snapshot):
        pass

    def _set_id(self, idx, obj):
        if self.reference_by_uuid:
            self.vars['uuid'][int(idx / 2)] = obj.__uuid__

    def _get_id(self, idx, obj):
        if self.reference_by_uuid:
            uuid = self.vars['uuid'][int(idx / 2)]
            if idx & 1:
                uuid = StorableObject.ruuid(uuid)

            obj.__uuid__ = uuid

    def load_indices(self):
        if self.reference_by_uuid:
            for pos, idx in enumerate(self.vars['index'][:]):
                self.index[idx] = pos

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

        self._set(idx, snapshot)

    def all(self):
        return peng.Trajectory(map(self.proxy, range(len(self))))

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
        st_idx = int(idx)
        self._set(st_idx, snapshot)

        return idx

    def idx(self, obj):
        return self.index[obj]

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
        return self.snapshot_class.__features__.classes

    @property
    def storables(self):
        return self.snapshot_class.__features__.storables

    def _set(self, idx, snapshot):
        [self.write(attr, idx, snapshot) for attr in self.storables]

    def _get(self, idx, snapshot):
        [setattr(snapshot, attr, self.vars[attr][idx]) for attr in self.storables]

    def initialize(self):
        super(FeatureSnapshotStore, self).initialize()

        for feature in self.classes:
            if hasattr(feature, 'netcdfplus_init'):
                feature.netcdfplus_init(self)


class FeatureSnapshotIndexedStore(BaseSnapshotIndexedStore):
    """
    An ObjectStore for Snapshots in netCDF files.
    """

    def __init__(self, descriptor):
        super(FeatureSnapshotIndexedStore, self).__init__(
            descriptor
    )

    @property
    def classes(self):
        return self.snapshot_class.__features__.classes

    @property
    def storables(self):
        return self.snapshot_class.__features__.storables

    def _set(self, idx, snapshot):
        [self.write(attr, idx, snapshot) for attr in self.storables]

    def _get(self, idx, snapshot):
        [setattr(snapshot, attr, self.vars[attr][idx]) for attr in self.storables]

    def initialize(self):
        super(FeatureSnapshotIndexedStore, self).initialize()

        for dim, size in self._dimensions.iteritems():
            self.storage.create_dimension(self.prefix + dim, size)

        for feature in self.classes:
            if hasattr(feature, 'netcdfplus_init'):
                feature.netcdfplus_init(self)

        self.storage.sync()


class SnapshotWrapperStore(ObjectStore):
    """
    A Store to store arbitrary snapshots
    """
    def __init__(self):
        super(SnapshotWrapperStore, self).__init__(peng.BaseSnapshot, json=False)

        self.type_list = {}
        self.store_list = []
        self.cv_list = {}

        self._treat_missing_snapshot_type = 'fail'


    @property
    def treat_missing_snapshot_type(self):
        return self._treat_missing_snapshot_type

    @treat_missing_snapshot_type.setter
    def treat_missing_snapshot_type(self, value):
        allowed = ['create', 'ignore', 'fail']
        if value not in allowed:
            raise ValueError('Only one of %s choices allowed.' % allowed)

        self._treat_missing_snapshot_type = value

    def _load(self, idx):
        store = self.vars['store'][idx]
        return store[idx]

    def _save(self, obj, idx):
        try:
            store, store_idx = self.type_list[obj.engine.descriptor]
            store[idx] = obj
            self.vars['store'][idx] = store_idx
            return store

        except KeyError:
            # Apparently there is no store yet to handle the given type of snapshot
            if self.treat_missing_snapshot_type == 'create':
                # we just create space for it
                store, store_idx = self.add_type(obj.engine.descriptor)
                store[idx] = obj
                self.vars['store'][idx] = store_idx
                return store

            elif self.treat_missing_snapshot_type == 'ignore':
                # we keep silent about it
                self.vars['store'][idx] = -1
                return None
            else:
                # we fail with cannot store
                raise RuntimeError(
                    (
                        'The store cannot hold snapshots of the given type : '
                        'class "%s" and dimensions %s. Try adding the snapshot type '
                        'using .add_type(snapshot).'
                    ) % (
                        obj.__class__.__name__,
                        obj.engine.descriptor.dimensions
                    )
                )

    def initialize(self):
        super(SnapshotWrapperStore, self).initialize()

        self.create_variable('store', 'index')

        self.storage.create_dimension('snapshottype')
        self.storage.create_dimension('cvcache')

        self.storage.create_variable('snapshottype', 'obj.stores', 'snapshottype')
        self.storage.create_variable('cvcache', 'obj.stores', 'cvcache')

    def add_type(self, descriptor):
        if isinstance(descriptor, peng.BaseSnapshot):
            descriptor = descriptor.engine.descriptor

        if descriptor in self.type_list:
            return self.type_list[descriptor]

        store = FeatureSnapshotIndexedStore(descriptor)

        store_idx = int(len(self.storage.dimensions['snapshottype']))
        store_name = 'snapshot' + str(store_idx)
        self.storage.register_store(store_name, store, False)

        # this will tell the store to add its own prefix for dimension names
        store.set_dimension_prefix_store(store)

        store.name = store_name
        self.storage.stores.save(store)

        self.type_list[descriptor] = (store, store_idx)
        self.store_list.append(store)
        self.storage.vars['snapshottype'][store_idx] = store

        self.storage.finalize_stores()
        self.storage.update_delegates()

        return store

    @staticmethod
    def _snapshot_store_name(idx):
        return 'snapshot' + str(idx)

    @staticmethod
    def to_descriptor(cls, dims):
        description = {}
        description.update(dims)
        description['class'] = cls

        return description

    def restore(self):
        for idx, store in enumerate(self.vars['snapshottype']):
            self.type_list[store.descriptor] = (store, idx)
            self.store_list.append(store)

    def save(self, obj, idx=None):
        # Cases to cover
        # A. uuid
        # 1. uuid not in self.index and ruuid not in self.index -> save
        # 2. uuid not in self.index and ruuid in self.index but not in store.index -> save
        # 3. uuid in self.index but not in store.index -> save
        # 4. uuid in self.index and in store.index -> return
        # 5. uuid not in self.index but ruuid in self.index and in store.index -> return
        # B. int
        # 1. obj not in self.index and obj._reversed not in self.index -> save
        # 2. obj not in self.index and obj._reversed in self.index but not in store.index -> save
        # 3. obj in self.index but not in store.index -> save
        # 4. obj in self.index and in store.index -> return
        # 5. obj not in self.index but obj._reversed in self.index and in store.index -> return
        if self.reference_by_uuid:

            if obj in self.index:
                # has could have been saved
                s_idx = self.index[obj]
                store_idx = self.variables['store'][s_idx]

                if not store_idx  == -1:
                    # let's see if the store has it stored
                    store = self.store_list[store_idx]
                    if s_idx in store.index:
                        # found it so return the correct reference
                        return self.reference(obj)
                    elif:

                else:
                    # if set to -1 then only an empty snapshot has been saved for CV purposes, etc
                    pass

            if hasattr(obj, '_idx'):
                if obj._store is self:
                    # is a proxy of a saved object so do nothing
                    return obj._idx
                else:
                    # it is stored but not in this store so we try storing the
                    # full snapshot which might be still in cache or memory
                    # if that is not the case it will be stored again. This can
                    # happen when you load from one store save to another. And load
                    # again after some time while the cache has been changed and try
                    # to save again the loaded object. We will not explicitly store
                    # a table that matches objects between different storages.
                    return self.save(obj.__subject__)
        else:
            if obj in self.index:
                # has could have been saved
                s_idx = self.index[obj]
                store_idx = self.variables['store'][s_idx]

                if not store_idx  == -1:
                    # let's see if the store has it stored
                    store = self.store_list[store_idx]
                    if s_idx in store.index:
                        # found it so return the correct reference
                        return self.reference(obj)
                    elif

                else:
                    # if set to -1 then only an empty snapshot has been saved for CV purposes, etc
                    pass

            if hasattr(obj, '_idx'):
                if obj._store is self:
                    # is a proxy of a saved object so do nothing
                    return obj._idx
                else:
                    # it is stored but not in this store so we try storing the
                    # full snapshot which might be still in cache or memory
                    # if that is not the case it will be stored again. This can
                    # happen when you load from one store save to another. And load
                    # again after some time while the cache has been changed and try
                    # to save again the loaded object. We will not explicitly store
                    # a table that matches objects between different storages.
                    return self.save(obj.__subject__)

        if not isinstance(obj, self.content_class):
            raise ValueError(
                'This store can only store object of base type "%s". Given obj is of type "%s". You'
                'might need to use another store.' % (self.content_class, obj.__class__.__name__)
            )

        n_idx = self.free()

        # mark as saved so circular dependencies will not result in infinite loops
        self.index[obj] = n_idx

        # make sure in nested saving that an IDX is not used twice!
        self.reserve_idx(n_idx)

        try:
            store = self._save(obj, n_idx)

            # store the name in the cache
            if hasattr(self, 'cache'):
                self.cache[n_idx] = obj

        except:
            # in case we did not succeed remove the mark as being saved
            del self.index[obj]
            self.release_idx(n_idx)
            raise

        self.release_idx(n_idx)
        self._set_id(n_idx, obj)

        return self.reference(obj)