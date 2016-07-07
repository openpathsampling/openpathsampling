import abc
from collections import OrderedDict

from openpathsampling.netcdfplus import StorableObject, LoaderProxy
from openpathsampling.netcdfplus.objects import UUIDDict, IndexedObjectStore
from openpathsampling.netcdfplus import NetCDFPlus, ObjectStore, LRUChunkLoadingCache
import openpathsampling.engines as peng


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


# =============================================================================================
# ABSTRACT BASE CLASS FOR SNAPSHOTS
# =============================================================================================

class BaseSnapshotStore(IndexedObjectStore):
    """
    An ObjectStore for Snapshots in netCDF files.
    """

    __metaclass__ = abc.ABCMeta

    def __init__(self, descriptor):
        """

        Attributes
        ----------
        descriptor : openpathsampling.engines.SnapshotDescriptor
            a descriptor knowing the snapshot class and a dictionary of
            dimensions and their lengths.

        """

        # Using a store with None as type will not interfere with the main snapshotstore
        super(BaseSnapshotStore, self).__init__(None, json=False)
        self.descriptor = descriptor
        self._dimensions = descriptor.dimensions
        self._cls = self.descriptor.snapshot_class

        # self._use_lazy_reversed = False
        # if hasattr(snapshot_class, '__features__'):
        #     if '_reversed' in snapshot_class.__features__.lazy:
        #         self._use_lazy_reversed = True

    def create_uuid_index(self):
        return UUIDReversalDict()

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
            self.prefix, self._cls.__name__, 'BaseSnapshot')

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
            'descriptor': self.descriptor,
        }

    def load(self, idx):
        pos = idx / 2

        # we want to load by uuid and it was not in cache.
        if pos in self.index:
            n_idx = self.index[pos]
        else:
            raise KeyError(idx)

        if n_idx < 0:
            return None

        # if it is in the cache, return it
        try:
            obj = self.cache[n_idx]
            if idx & 1:
                obj = obj.reversed

            return obj

        except KeyError:
            pass

        obj = self._load(n_idx)
        self.cache[n_idx] = obj

        if idx & 1:
            obj = obj.reversed

        return obj

    def _load(self, idx):
        """
        Load a snapshot from the storage.

        Parameters
        ----------
        idx : int
            the integer index of the snapshot to be loaded

        Returns
        -------
        snapshot : :obj:`openpathsampling.engines.BaseSnapshot`
            the loaded snapshot instance
        """

        # if not load and return it
        st_idx = int(idx)

        obj = self._cls.__new__(self._cls)
        self._cls.init_empty(obj)
        self._get(st_idx, obj)
        return obj

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

    @abc.abstractmethod
    def _get(self, idx, snapshot):
        pass

    @abc.abstractmethod
    def _set(self, idx, snapshot):
        pass

    def _get_id(self, idx, obj):
        if self.reference_by_uuid:
            uuid = self.vars['uuid'][int(idx / 2)]
            if idx & 1:
                uuid = StorableObject.ruuid(uuid)

            obj.__uuid__ = uuid

    def _set_id(self, idx, obj):
        if self.reference_by_uuid:
            self.vars['uuid'][int(idx / 2)] = obj.__uuid__

    def load_indices(self):
        if self.reference_by_uuid:
            for pos, idx in enumerate(self.vars['index'][:]):
                self.index[idx] = pos

    def all(self):
        return peng.Trajectory(map(self.proxy, range(len(self))))

    def duplicate(self, snapshot):
        """
        Store a duplicate of the snapshot as new

        Parameters
        ----------
        snapshot : :class:`openpathsampling.engines.BaseSnapshot`

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

    def __init__(self, descriptor):
        super(FeatureSnapshotStore, self).__init__(descriptor)

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
        self.store_snapshot_list = []
        self.store_cv_list = []
        self.cv_list = {}

        # default way to handle unknown snapshot types is to create
        # a single store for the first type tried to be stored
        # if you want to store more than one snapshots you
        # need to add them manually
        self._treat_missing_snapshot_type = 'single'

        # if set to true snapshots will not be stored but merely registered
        # so CVs will be storable
        self.only_mention = False

    @property
    def treat_missing_snapshot_type(self):
        return self._treat_missing_snapshot_type

    @treat_missing_snapshot_type.setter
    def treat_missing_snapshot_type(self, value):
        allowed = ['create', 'ignore', 'fail', 'single']
        if value not in allowed:
            raise ValueError('Only one of %s choices allowed.' % allowed)

        self._treat_missing_snapshot_type = value

    def _load(self, idx):
        store_idx = self.vars['store'][idx / 2]

        if store_idx < 0:
            raise KeyError('IDX "' + idx + '" not found in storage')
        else:
            store = self.store_snapshot_list[store_idx]
            return store[idx]

    def __len__(self):
        return len(self.storage.dimensions[self.prefix]) * 2

    def initialize(self):
        super(SnapshotWrapperStore, self).initialize()

        self.create_variable('store', 'index')

        self.storage.create_dimension('snapshottype')
        self.storage.create_dimension('cvcache')

        self.storage.create_variable('snapshottype', 'obj.stores', 'snapshottype')
        self.storage.create_variable('cvcache', 'obj.stores', 'cvcache')

    def add_type(self, descriptor):
        if isinstance(descriptor, peng.BaseSnapshot):
            template = descriptor
            descriptor = descriptor.engine.descriptor
        else:
            template = None

        if descriptor in self.type_list:
            return self.type_list[descriptor]

        store = FeatureSnapshotStore(descriptor)

        store_idx = int(len(self.storage.dimensions['snapshottype']))
        store_name = 'snapshot' + str(store_idx)
        self.storage.register_store(store_name, store, False)

        # this will tell the store to add its own prefix for dimension names
        store.set_dimension_prefix_store(store)

        store.name = store_name
        self.storage.stores.save(store)

        self.type_list[descriptor] = (store, store_idx)
        self.store_snapshot_list.append(store)
        self.storage.vars['snapshottype'][store_idx] = store

        self.storage.finalize_stores()
        self.storage.update_delegates()

        if template:
            self.save(template)

        return store, store_idx

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
        for idx, store in enumerate(self.storage.vars['snapshottype']):
            self.type_list[store.descriptor] = (store, idx)
            self.store_snapshot_list.append(store)

        for idx, store in enumerate(self.storage.vars['cvcache']):
            cv = self.storage.cvs._load(int(store.name[2:]))
            self.cv_list[cv] = (store, idx)

        self.load_indices()

    def load_indices(self):
        if self.reference_by_uuid:
            for idx, uuid in enumerate(self.vars['uuid'][:]):
                self.index[uuid] = idx * 2

    def get_cv_cache(self, idx):
        store_name = SnapshotWrapperStore._get_cv_name(idx)
        if store_name in self.storage.stores.name_idx:
            return self.storage.stores[store_name]
        else:
            return None

    def save(self, obj, idx=None):
        n_idx = None

        if self.reference_by_uuid:
            if obj in self.index:
                n_idx = self.index[obj]
        else:
            if hasattr(obj, '_idx'):
                if obj._store is self:
                    # is a proxy of a saved object so do nothing
                    return obj._idx
                else:
                    return self.save(obj.__subject__)

            if obj in self.index:
                n_idx = self.index[obj]

            elif obj._reversed is not None:
                # if the object has no reversed present, then the reversed does not
                # exist yet and hence it cannot be in the index, so no checking
                if obj._reversed in self.index:
                    n_idx = self.index[obj._reversed] ^ 1

        store_idx = None

        if n_idx is not None:
            # snapshot is mentioned
            store_idx = self.variables['store'][n_idx / 2]
            if not store_idx == -1:
                # and stored
                if self.reference_by_uuid:
                    return self.reference(obj)
                else:
                    return n_idx

        else:
            if self.only_mention:
                # only mention but not really store snapshots
                self.vars['store'][n_idx / 2] = -1
                n_idx = self.free()
                self._set_id(n_idx, obj)
                return self.reference(obj)

        if not isinstance(obj, self.content_class):
            raise ValueError(
                'This store can only store object of base type "%s". Given obj is of type "%s". You'
                'might need to use another store.' % (self.content_class, obj.__class__.__name__)
            )

        if n_idx is None:
            n_idx = self.free()

            self._save(obj, n_idx)
            self._complete_cv(obj, n_idx)
            self._set_id(n_idx, obj)
        else:
            self.vars['store'][n_idx / 2] = store_idx
            self.index[obj] = n_idx
            store = self.store_snapshot_list[store_idx]
            store[n_idx / 2] = obj

        self.cache[n_idx] = obj

        return self.reference(obj)

    def _complete_cv(self, obj, pos):
        for cv, (cv_store, cv_idx) in self.cv_list.items():
            if cv_store.auto_complete:
                value = cv._cache_dict._get(obj)
                if value is None:
                    # not in cache so compute it if possible
                    if cv._eval_dict:
                        value = cv._eval_dict([obj])[0]
                    else:
                        value = None

                if value is not None:
                    if cv_store.allow_partial:
                        cv_store[obj] = value
                    else:
                        if cv_store.time_reversible:
                            n_idx = pos / 2
                        else:
                            n_idx = pos

                        cv_store.vars['value'][n_idx] = value
                        cv_store.cache[n_idx] = value

    def _save(self, obj, idx):
        try:
            store, store_idx = self.type_list[obj.engine.descriptor]
            self.vars['store'][idx / 2] = store_idx
            self.index[obj] = idx
            store[idx / 2] = obj
            return store

        except KeyError:
            # Apparently there is no store yet to handle the given type of snapshot
            mode = self.treat_missing_snapshot_type
            if mode == 'create' or \
                    (mode == 'single' and len(self.storage.dimensions['snapshottype']) == 0):
                # we just create space for it
                store, store_idx = self.add_type(obj.engine.descriptor)
                self.vars['store'][idx / 2] = store_idx
                self.index[obj] = idx
                store[idx / 2] = obj
                return store

            elif self.treat_missing_snapshot_type == 'ignore':
                # we keep silent about it
                self.vars['store'][idx / 2] = -1
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

    def free(self):
        idx = len(self)
        while idx in self._free:
            idx += 2

        return idx

    def get_uuid_index(self, obj):
        n_idx = None

        if self.reference_by_uuid:
            if obj in self.index:
                n_idx = self.index[obj]
        else:
            if hasattr(obj, '_idx'):
                if obj._store is self:
                    # is a proxy of a saved object so do nothing
                    return obj._idx

            if obj in self.index:
                n_idx = self.index[obj]

            elif obj._reversed is not None:
                if obj._reversed in self.index:
                    n_idx = self.index[obj._reversed] ^ 1

        if n_idx is None:
            # if the obj is not know, add it to the file and index, but
            # store only a reference and not the full object
            # this can later be done using .save(obj)
            n_idx = self.free()
            self.variables['store'][n_idx / 2] = -1
            self.index[obj] = n_idx
            self._set_id(n_idx, obj)

    @staticmethod
    def _get_cv_name(cv_idx):
        return 'cv' + str(cv_idx)

    def add_cv(self, cv, template, allow_partial=None, auto_complete=None, chunksize=None):
        """

        Parameters
        ----------
        cv : :obj:`openpathsampling.CollectiveVariable`
        template : :obj:`openpathsampling.engines.BaseSnapshot`
        chunksize : int
        allow_partial : bool
        auto_complete : bool

        Returns
        -------
        :obj:`openpathsampling.netcdfplus.ObjectStore`
        int
        """
        if cv in self.cv_list:
            return self.cv_list[cv]

        if allow_partial is None:
            allow_partial = cv.diskcache_allow_partial
        if auto_complete is None:
            auto_complete = cv.diskcache_auto_complete
        if chunksize is None:
            chunksize = cv.diskcache_chunksize
        if template is None:
            template = cv.diskcache_template

        time_reversible = cv.cv_time_reversible

        if not time_reversible:
            # in the rare case of not time_reversible we need to store partial for now
            allow_partial = True

        if not allow_partial:
            # in complete mode we need to force chunksize one to match it to snapshots
            chunksize = 1

            # and autocomplete is on in complete
            auto_complete = True

        # determine value type and shape
        params = NetCDFPlus.get_value_parameters(cv(template))
        shape = params['dimensions']

        if shape is None:
            chunksizes = None
        else:
            chunksizes = tuple(params['dimensions'])

        cv_idx = self.storage.cvs.index[cv]
        store = SnapshotValueStore(
            time_reversible=time_reversible,
            allow_partial=allow_partial,
            auto_complete=auto_complete,
            chunksize=chunksize
        )

        store_name = SnapshotWrapperStore._get_cv_name(cv_idx)
        self.storage.create_store(store_name, store, False)

        if store.allow_partial:
            # we are not using the .initialize function here since we
            # only have one variable and only here know its shape
            self.storage.create_dimension(store.prefix, 0)

            if shape is not None:
                shape = tuple(list(shape))
                chunksizes = tuple([chunksize] + list(chunksizes))
            else:
                shape = tuple()
                chunksizes = tuple([chunksize])

            # create the variable
            store.create_variable(
                'value',
                var_type=params['var_type'],
                dimensions=shape,
                chunksizes=chunksizes,
                simtk_unit=params['simtk_unit'],
            )

            store.create_variable('index', 'index')

        else:
            chunksize = 1
            if shape is not None:
                shape = tuple(['snapshots'] + list(shape))
                chunksizes = tuple([chunksize] + list(chunksizes))
            else:
                shape = tuple(['snapshots'])
                chunksizes = tuple([chunksize])

            # create the variable
            store.storage.create_variable(
                store_name + '_value',
                var_type=params['var_type'],
                dimensions=shape,
                chunksizes=chunksizes,
                simtk_unit=params['simtk_unit'],
            )

        setattr(store, 'value', self.storage.vars[store_name + '_value'])

        store_idx = int(len(self.storage.dimensions['cvcache']))
        self.cv_list[cv] = (store, store_idx)
        self.storage.vars['cvcache'][store_idx] = store

        # use the cache and function of the CV to fill the store when it is made
        if not allow_partial:

            if self.reference_by_uuid:
                indices = self.vars['uuid'][:]
            else:
                indices = range(0, len(self), 2)

            for pos, idx in enumerate(indices):

                proxy = LoaderProxy(self.storage.snapshots, idx)
                value = cv._cache_dict._get(proxy)

                if value is None:
                    # not in cache so compute it if possible
                    if cv._eval_dict:
                        value = cv._eval_dict([proxy])[0]
                    else:
                        value = None

                if value is not None:
                    if time_reversible:
                        n_idx = pos / 2
                    else:
                        n_idx = pos

                    store.vars['value'][n_idx] = value
                    store.cache[n_idx] = value

        cv.set_cache_store(store)
        return store, store_idx

    def create_uuid_index(self):
        return UUIDReversalDict()

    def _get_id(self, idx, obj):
        if self.reference_by_uuid:
            uuid = self.vars['uuid'][int(idx / 2)]
            if idx & 1:
                if obj._reversed:
                    obj._reversed.__uuid__ = uuid

                uuid = StorableObject.ruuid(uuid)

            obj.__uuid__ = uuid

    def _set_id(self, idx, obj):
        if self.reference_by_uuid:
            self.vars['uuid'][idx / 2] = obj.__uuid__

    def complete_cv(self, cv):
        """
        Will fill in all missing values for a CV and store these.

        Parameters
        ----------
        cv

        Returns
        -------

        """
        if self.reference_by_uuid:
            if 'partial':
                uuids = self.vars['uuid'][:]
                stores = self.vars['store'][:]
                cv_store, cv_index = self.cv_list[cv]
                for pos, store in enumerate(stores):
                    uuid = uuids[pos]
                    if store >= 0:
                        cv_store[uuid] = cv(self.load(uuid))
            elif 'complete':
                pass
        else:
            stores = self.vars['store'][:]
            cv_store, cv_index = self.cv_list[cv]
            for pos, store in enumerate(stores):
                if store >= 0:
                    cv_store[pos] = cv(self[pos])

    def idx(self, obj):
        """
        Return the index in this store for a given object

        Parameters
        ----------
        obj : :py:class:`openpathsampling.engines.Snapshot`
            the object that can be stored in this store for which its index is
            to be returned

        Returns
        -------
        int or None
            The integer index of the given object or None if it is not stored yet
        """
        try:
            return self.index[obj]
        except KeyError:
            try:
                return self.index[obj._reversed] ^ 1
            except KeyError:
                raise KeyError(obj)

    def cache_all_cvs(self):
        for store, idx in self.cv_list.values():
            store.fill_cache()

    def pos(self, obj):
        if hasattr(obj, '_idx'):
            if self.reference_by_uuid:
                pos = self.index.get(obj)
            else:
                pos = obj._idx
        else:
            pos = self.index.get(obj)

            if pos is None and not self.reference_by_uuid:
                if obj._reversed:
                    pos = self.index.get(obj._reversed)
                    if pos is None:
                        return None

                    pos ^= 1
                else:
                    return None

        return pos


class SnapshotValueStore(ObjectStore):
    def __init__(
            self,
            time_reversible=True,
            allow_partial=False,
            auto_complete=True,
            chunksize=100
    ):
        super(SnapshotValueStore, self).__init__(None)
        self.snapshot_index = None
        if not time_reversible and not allow_partial:
            raise RuntimeError(
                'Only time_reversible CVs can currently be stored using mode "complete"')

        self.time_reversible = time_reversible
        self.allow_partial = allow_partial
        self.auto_complete = auto_complete
        self.chunksize = chunksize

    def to_dict(self):
        return {
            'time_reversible': self.time_reversible,
            'allow_partial': self.allow_partial,
            'auto_complete': self.auto_complete
        }

    def create_uuid_index(self):
        return dict()

    def create_int_index(self):
        return dict()

    def register(self, storage, prefix):
        super(SnapshotValueStore, self).register(storage, prefix)
        self.snapshot_pos = self.storage.snapshots.pos

    def __len__(self):
        return len(self.variables['value'])

    # =============================================================================
    # LOAD/SAVE DECORATORS FOR CACHE HANDLING
    # =============================================================================

    def load(self, idx):
        pos = self.snapshot_pos(idx)

        if pos is None:
            if self.reference_by_uuid:
                return None
            else:
                if idx._reversed:
                    pos = self.snapshot_pos(idx._reversed)

                    if pos is None:
                        return None

                    pos ^= 1
                else:
                    return None

        if self.time_reversible:
            pos /= 2

        if self.allow_partial:
            # we want to load by uuid and it was not in cache.
            if pos in self.index:
                n_idx = self.index[pos]
            else:
                return None
            if n_idx < 0:
                return None
        else:
            if pos < len(self):
                n_idx = pos
            else:
                return None

        # if it is in the cache, return it
        try:
            obj = self.cache[n_idx]
            return obj

        except KeyError:
            pass

        obj = self.vars['value'][n_idx]

        self.cache[n_idx] = obj

        return obj

    def __setitem__(self, idx, value):
        pos = self.snapshot_pos(idx)

        if pos is None:
            return

        if self.time_reversible:
            pos /= 2

        if self.allow_partial:
            if pos in self.index:
                return

            n_idx = self.free()
        else:
            if pos < len(self):
                return

            n_idx = idx

        if self.allow_partial:
            # only if partial storage is used store index and update
            self.vars['index'][n_idx] = pos
            self.index[pos] = n_idx

        self.vars['value'][n_idx] = value
        self.cache[n_idx] = value

    def fill_cache(self):
        self.cache.load_max()

    def restore(self):
        if self.allow_partial:  # only if partial storage is used
            for pos, idx in enumerate(self.vars['index'][:]):
                self.index[idx] = pos

        self.set_caching(LRUChunkLoadingCache(
            chunksize=self.chunksize,
            max_chunks=1000,
            variable=self.vars['value']
        ))

    def initialize(self):
        pass

    def __getitem__(self, item):
        # enable numpy style selection of objects in the store
        try:
            if isinstance(item, peng.BaseSnapshot):
                return self.load(item)
            elif type(item) is list:
                return [self.load(idx) for idx in item]
            elif item is Ellipsis:
                return self.iterator()
        except KeyError:
            return None

    def get(self, item):
        if self.allow_partial:
            if item in self.index:
                return self[item]
            else:
                return None
        else:
            return self[item]
