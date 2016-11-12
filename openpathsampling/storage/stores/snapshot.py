import abc
from collections import OrderedDict

from openpathsampling.netcdfplus import StorableObject, LoaderProxy
from openpathsampling.netcdfplus.objects import UUIDDict, IndexedObjectStore
from openpathsampling.netcdfplus import NetCDFPlus, ObjectStore, \
    LRUChunkLoadingCache
import openpathsampling.engines as peng

from openpathsampling.netcdfplus import with_timing_logging

import logging
from uuid import UUID

logger = logging.getLogger(__name__)
init_log = logging.getLogger('openpathsampling.initialization')


class UUIDReversalDict(UUIDDict):
    @staticmethod
    def rev_id(obj):
        return StorableObject.ruuid(UUIDReversalDict.id(obj))

    def __setitem__(self, key, value, **kwargs):
        OrderedDict.__setitem__(self, self.id(key), value)
        OrderedDict.__setitem__(self, self.rev_id(key), value ^ 1)

    def __delitem__(self, key, **kwargs):
        OrderedDict.__delitem__(self, self.id(key))
        OrderedDict.__delitem__(self, self.rev_id(key))


# ==============================================================================
# ABSTRACT BASE CLASS FOR SNAPSHOTS
# ==============================================================================

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

        # Using a store with None as type will not interfere with the main
        # `SnapshotWrapperStore`
        super(BaseSnapshotStore, self).__init__(None, json=False)
        self.descriptor = descriptor
        self._dimensions = descriptor.dimensions
        self._cls = self.descriptor.snapshot_class

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
        This make storing CVs easier. This function allows to get the paired
        index or the index of snapshot.reversed

        The implementation uses the trick that all you have to do is flip the
        lowest bit that determines even or odd.

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

        # # if it is in the cache, return it
        # try:
        #     obj = self.cache[n_idx]
        #     if idx & 1:
        #         obj = obj.reversed
        #
        #     return obj
        #
        # except KeyError:
        #     pass

        obj = self._load(n_idx)

        # self.cache[n_idx] = obj

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

    def save(self, obj, idx=None):
        pos = idx / 2

        if pos in self.index:
            # has been saved so quit and do nothing
            return idx

        # n_idx = self.free() / 2
        n_idx = len(self.index) / 2

        # mark as saved so circular dependencies will not cause infinite loops
        self.index.append(pos)

        # make sure in nested saving that an IDX is not used twice!
        # self.reserve_idx(n_idx)

        logger.debug('Saving ' + str(type(obj)) + ' using IDX #' + str(n_idx))

        try:
            self._save(obj, n_idx)
            self.vars['index'][n_idx] = pos

            # store the name in the cache
            # if hasattr(self, 'cache'):
            self.cache[n_idx] = obj

        except:
            logger.debug('Problem saving %d !' % n_idx)
            # in case we did not succeed remove the mark as being saved
            del self.index[pos]
            # self.release_idx(n_idx)
            raise

        # self.release_idx(n_idx)
        self._set_id(n_idx, obj)

        return idx

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
        uuid = self.vars['uuid'][int(idx / 2)]

        if idx & 1:
            uuid = StorableObject.ruuid(uuid)

        obj.__uuid__ = uuid

    def _set_id(self, idx, obj):
        self.vars['uuid'][int(idx / 2)] = obj.__uuid__

    def load_indices(self):
        self.index.extend(self.vars['index'])

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
            the index used for storing it in the store. This is the
            save as used by save.

        Notes
        -----
        This will circumvent the caching and indexing completely. This would be
        equivalent of creating a copy of the current snapshot and store this one
        and throw the copy away, leaving the given snapshot untouched. This
        allows you to treat the snapshot as mutual.

        The use becomes more obvious when applying to storing trajectories.
        The only way to make use of this feature is using the returned `idx`

        >>> idx = store.duplicate(snap)
        >>> loaded = store[idx]  # return a duplicated as new object
        >>> proxy = paths.LoaderProxy(store, idx) # use duplicate w/o loading

        """
        idx = self.free()
        st_idx = int(idx)
        self._set(st_idx, snapshot)

        return idx

    def __iter__(self):
        for idx in range(0, len(self), 2):
            yield self[idx]

    def __getitem__(self, item):
        """
        Enable numpy style selection of object in the store
        """
        try:
            if type(item) is int or type(item) is str or type(item) is UUID:
                return self.load(item)
            elif type(item) is slice:
                return [self.load(idx)
                        for idx in range(*item.indices(len(self)))]
            elif type(item) is list:
                return [self.load(idx) for idx in item]
            elif item is Ellipsis:
                return iter(self)
        except KeyError:
            return None

    def __len__(self):
        return len(self.storage.dimensions[self.prefix]) * 2


# ==============================================================================
# FEATURE BASED SINGLE CLASS FOR ALL SNAPSHOT TYPES
# ==============================================================================

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
        [setattr(snapshot, attr, self.vars[attr][idx])
         for attr in self.storables]

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
        super(SnapshotWrapperStore, self).__init__(
            peng.BaseSnapshot,
            json=False
        )

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

    def load(self, idx):
        """
        Returns an object from the storage.

        Parameters
        ----------
        idx : int
            the integer index of the object to be loaded

        Returns
        -------
        :py:class:`openpathsampling.netcdfplus.base.StorableObject`
            the loaded object
        """

        if type(idx) is long:
            # we want to load by uuid and it was not in cache.
            print 'long', idx
            if idx in self.index:
                n_idx = self.index[idx]
            else:
                if self.fallback_store is not None:
                    return self.fallback_store.load(idx)
                elif self.storage.fallback is not None:
                    return self.storage.fallback.stores[self.name].load(idx)
                else:
                    raise ValueError(
                        'str %s not found in storage or fallback' % idx)

        elif type(idx) is not int:
            raise ValueError(
                ('indices of type "%s" are not allowed in named storage '
                 '(only str and int)') % type(idx).__name__
            )
        else:
            n_idx = int(idx)

        print n_idx

        if n_idx < 0:
            return None

        # if it is in the cache, return it
        try:
            obj = self.cache[n_idx]
            logger.debug('Found IDX #' + str(idx) + ' in cache. Not loading!')
            return obj

        except KeyError:
            try:
                obj = self.cache[n_idx ^ 1].reversed
                logger.debug('Found IDX #' + str(idx) +
                             ' reversed in cache. Not loading!')
                return obj
            except KeyError:
                pass

        logger.debug(
            'Calling load object of type ' + self.content_class.__name__ +
            ' and IDX #' + str(idx))

        if n_idx >= len(self):
            logger.warning(
                'Trying to load from IDX #' + str(n_idx) +
                ' > number of object ' + str(len(self)))
            return None
        elif n_idx < 0:
            logger.warning(
                'Trying to load negative IDX #' + str(n_idx) +
                ' < 0. This should never happen!!!')
            raise RuntimeError(
                'Loading of negative int should result in no object. '
                'This should never happen!')
        else:
            obj = self._load(n_idx)

        logger.debug(
            'Calling load object of type %s and IDX # %d ... DONE' %
            (self.content_class.__name__, n_idx))

        if obj is not None:
            self._get_id(n_idx, obj)
            # self.index[obj.__uuid__] = n_idx
            self.cache[n_idx] = obj

            logger.debug(
                'Try loading UUID object of type %s and IDX # %d ... DONE' %
                (self.content_class.__name__, n_idx))

        logger.debug(
            'Finished load object of type %s and IDX # %d ... DONE' %
            (self.content_class.__name__, n_idx))

        return obj

    def _load(self, idx):
        store_idx = int(self.variables['store'][idx / 2])

        if store_idx < 0:
            if self.fallback_store is not None:
                return self.fallback_store.load(idx)
            elif self.storage.fallback is not None:
                return self.storage.fallback.snapshots.load(idx)
            else:
                raise KeyError(
                    'str %s not found in storage or fallback' % idx)
        else:
            store = self.store_snapshot_list[store_idx]
            snap = store[idx]
            return snap

    def __len__(self):
        return len(self.storage.dimensions[self.prefix]) * 2

    def initialize(self):
        super(SnapshotWrapperStore, self).initialize()

        self.create_variable('store', 'index')

        self.storage.create_dimension('snapshottype')
        self.storage.create_dimension('cvcache')

        self.storage.create_variable(
            'snapshottype',
            'obj.stores',
            'snapshottype')
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

        self.storage.cvs.load_indices()

        for idx, store in enumerate(self.storage.vars['cvcache']):
            cv_st_idx = int(store.name[2:])

            cv = self.storage.cvs[self.storage.cvs.vars['uuid'][cv_st_idx]]
            self.cv_list[cv] = (store, idx)

        self.load_indices()

    @with_timing_logging
    def load_indices(self):
        # TODO: Update with ReversedHashedList
        for idx, uuid in enumerate(self.vars['uuid'][:]):
            self.index[uuid] = idx * 2

    def get_cv_cache(self, idx):
        store_name = SnapshotWrapperStore._get_cv_name(idx)

        if store_name in self.storage.stores.name_idx:
            store = self.storage.stores[store_name]
            return store
        else:
            return None

    def mention(self, snapshot):
        """
        Save a shallow copy

        Parameters
        ----------
        snapshot : `openpathsampling.engines.BaseSnapshot`
            the snapshot to be shallow stored

        Returns
        -------
        int or `UUID`
            the reference to find the shallow! object

        """
        current_mention = self.only_mention
        self.only_mention = True
        ref = self.save(snapshot)
        self.only_mention = current_mention
        return ref

    def save(self, obj, idx=None):
        n_idx = None

        if obj.__uuid__ in self.index:
            n_idx = self.index[obj.__uuid__]

        if n_idx is not None:
            # snapshot is mentioned
            store_idx = int(self.variables['store'][n_idx / 2])
            if not store_idx == -1:
                # and stored
                return self.reference(obj)

        if self.only_mention:
            if n_idx is None:
                n_idx = self.free()

                # only mention but not really store snapshots
                self.vars['store'][n_idx / 2] = -1
                self.index[obj.__uuid__] = n_idx
                self._auto_complete_single_snapshot(obj, n_idx)
                self._set_id(n_idx, obj)

            return self.reference(obj)

        if not isinstance(obj, self.content_class):
            raise ValueError(
                ('This store can only store object of base type "%s". '
                 'Given obj is of type "%s". You'
                 'might need to use another store.') %
                (self.content_class, obj.__class__.__name__)
            )

        if n_idx is None:
            n_idx = self.free()

            self._save(obj, n_idx)
            self._auto_complete_single_snapshot(obj, n_idx)
            self._set_id(n_idx, obj)
        else:
            store, store_idx = self.type_list[obj.engine.descriptor]
            self.vars['store'][n_idx / 2] = store_idx
            store[n_idx] = obj

        self.cache[n_idx] = obj

        return self.reference(obj)

    def _save(self, obj, n_idx):
        try:
            store, store_idx = self.type_list[obj.engine.descriptor]
            self.vars['store'][n_idx / 2] = store_idx
            self.index[obj.__uuid__] = n_idx
            store[n_idx] = obj
            return store

        except KeyError:
            # there is no store yet to handle the given type of snapshot
            mode = self.treat_missing_snapshot_type
            if mode == 'create' or \
                    (mode == 'single' and
                             len(self.storage.dimensions['snapshottype']) == 0):
                # we just create space for it
                store, store_idx = self.add_type(obj.engine.descriptor)
                self.vars['store'][n_idx / 2] = store_idx
                self.index[obj.__uuid__] = n_idx
                store[n_idx] = obj
                return store

            elif self.treat_missing_snapshot_type == 'ignore':
                # we keep silent about it
                self.vars['store'][n_idx / 2] = -1
                return None
            else:
                # we fail with cannot store
                raise RuntimeError(
                    (
                        'The store cannot hold snapshots of the given type : '
                        'class "%s" and dimensions %s. Try adding the '
                        'snapshot type using .add_type(snapshot).'
                    ) % (
                        obj.__class__.__name__,
                        obj.engine.descriptor.dimensions
                    )
                )

    def _auto_complete_single_snapshot(self, obj, pos):
        for cv, (cv_store, cv_idx) in self.cv_list.items():
            if not cv_store.allow_incomplete:
                value = cv._cache_dict._get(obj)
                if value is None:
                    # not in cache so compute it if possible
                    if cv._eval_dict:
                        value = cv._eval_dict([obj])[0]

                if value is not None:
                    if cv_store.allow_incomplete:
                        cv_store[obj] = value
                    else:
                        if cv_store.time_reversible:
                            n_idx = pos / 2
                        else:
                            n_idx = pos

                        cv_store.vars['value'][n_idx] = value
                        cv_store.cache[n_idx] = value

    def complete_cv(self, cv):
        """
        Compute all missing values of a CV and store them


        Parameters
        ----------
        cv : :obj:`openpathsampling.CollectiveVariable`


        """
        if cv not in self.cv_list:
            return

        cv_store = self.cv_list[cv][0]

        if cv_store.allow_incomplete:
            # for complete this does not make sense

            # TODO: Make better looping over this to not have
            # to load all the indices at once
            # can be problematic for 10M+ stored snapshots
            indices = self.vars['uuid'][:]

            for pos, idx in enumerate(indices):
                if not cv_store.time_reversible:
                    pos *= 2

                proxy = None

                if pos not in cv_store.index:
                    # this value is not stored to go ahead

                    proxy = self.storage.snapshots[idx]

                    # get from cache first, this is fastest
                    value = cv._cache_dict._get(proxy)

                    if value is None:
                        # not in cache so compute it if possible
                        if cv._eval_dict:
                            value = cv._eval_dict([proxy])[0]
                        else:
                            value = None

                    if value is not None:
                        n_idx = cv_store.free()

                        cv_store.vars['value'][n_idx] = value
                        cv_store.vars['index'][n_idx] = pos
                        cv_store.index[pos] = n_idx
                        cv_store.cache[n_idx] = value

                if not cv_store.time_reversible:
                    pos += 1
                    if pos not in cv_store.index:
                        if proxy is None:
                            proxy = self.storage.snapshots[idx]

                        if proxy._reversed is not None:
                            proxy = proxy._reversed
                        else:
                            proxy = proxy.reversed

                        # get from cache first, this is fastest
                        value = cv._cache_dict._get(proxy)

                        if value is None:
                            # not in cache so compute it if possible
                            if cv._eval_dict:
                                value = cv._eval_dict([proxy])[0]
                            else:
                                value = None

                        if value is not None:
                            n_idx = cv_store.free()

                            cv_store.vars['value'][n_idx] = value
                            cv_store.vars['index'][n_idx] = pos
                            cv_store.index[pos] = n_idx
                            cv_store.cache[n_idx] = value

    def sync_cv(self, cv):
        """
        Store all cached values of a CV in the diskcache

        Parameters
        ----------
        cv : :obj:`openpathsampling.CollectiveVariable`


        """

        if cv not in self.cv_list:
            return

        cv_store = self.cv_list[cv][0]

        # for complete this does not make sense
        if cv_store.allow_incomplete:

            # loop all objects in the fast CV cache
            for obj, value in cv._cache_dict.cache.iteritems():
                if value is not None:
                    pos = self.pos(obj)

                    # if the snapshot is not saved, nothing we can do
                    if pos is None:
                        continue

                    if cv_store.time_reversible:
                        pos /= 2

                    if pos in cv_store.index:
                        # this value is stored so skip it
                        continue

                    n_idx = cv_store.free()

                    cv_store.vars['value'][n_idx] = value
                    cv_store.vars['index'][n_idx] = pos
                    cv_store.index[pos] = n_idx
                    cv_store.cache[n_idx] = value

    def free(self):
        idx = len(self)
        while idx in self._free:
            idx += 2

        return idx

    def get_uuid_index(self, obj):
        n_idx = None

        if obj.__uuid__ in self.index:
            n_idx = self.index[obj.__uuid__]

        if n_idx is None:
            # if the obj is not know, add it to the file and index, but
            # store only a reference and not the full object
            # this can later be done using .save(obj)
            n_idx = self.free()
            self.variables['store'][n_idx / 2] = -1
            self.index[obj.__uuid__] = n_idx
            self._set_id(n_idx, obj)

    @staticmethod
    def _get_cv_name(cv_idx):
        return 'cv' + str(cv_idx)

    def add_cv(self, cv, template, allow_incomplete=None, chunksize=None):
        """

        Parameters
        ----------
        cv : :obj:`openpathsampling.CollectiveVariable`
        template : :obj:`openpathsampling.engines.BaseSnapshot`
        chunksize : int
        allow_incomplete : bool

        Returns
        -------
        :obj:`openpathsampling.netcdfplus.ObjectStore`
        int
        """
        if cv in self.cv_list:
            return self.cv_list[cv]

        if allow_incomplete is None:
            allow_incomplete = cv.diskcache_allow_incomplete
        if chunksize is None:
            chunksize = cv.diskcache_chunksize
        if template is None:
            template = cv.diskcache_template

        time_reversible = cv.cv_time_reversible

        if not time_reversible:
            # in the rare case of not time_reversible we use store partial
            allow_incomplete = True

        if not allow_incomplete:
            # in complete mode we force chunk size one to match it to snapshots
            chunksize = self.default_store_chunk_size

        # determine value type and shape
        params = NetCDFPlus.get_value_parameters(cv(template))
        shape = params['dimensions']

        if shape is None:
            chunksizes = None
        else:
            chunksizes = tuple(params['dimensions'])

        cv_idx = self.storage.cvs.index[cv.__uuid__]
        store = SnapshotValueStore(
            time_reversible=time_reversible,
            allow_incomplete=allow_incomplete,
            chunksize=chunksize
        )

        store_name = SnapshotWrapperStore._get_cv_name(cv_idx)
        self.storage.create_store(store_name, store, False)

        if store.allow_incomplete:
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
            chunksize = self.default_store_chunk_size
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

        store.initialize()

        store_idx = int(len(self.storage.dimensions['cvcache']))
        self.cv_list[cv] = (store, store_idx)
        self.storage.vars['cvcache'][store_idx] = store

        # use the cache and function of the CV to fill the store when it is made
        if not allow_incomplete:

            indices = self.vars['uuid'][:]

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
                    store.vars['value'][pos] = value
                    store.cache[pos] = value

        cv.set_cache_store(store)
        return store, store_idx

    def create_uuid_index(self):
        return UUIDReversalDict()

    def _get_id(self, idx, obj):
        uuid = self.vars['uuid'][int(idx / 2)]

        if idx & 1:
            uuid = StorableObject.ruuid(uuid)

        obj.__uuid__ = uuid
        if idx & 1:
            if obj._reversed:
                obj._reversed.__uuid__ = uuid

            uuid = StorableObject.ruuid(uuid)
            obj.__uuid__ = uuid

    def _set_id(self, idx, obj):
        self.vars['uuid'][idx / 2] = obj.__uuid__

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
            The integer index of the given object or None if it is not
            stored yet
        """
        try:
            return self.index[obj.__uuid__]
        except KeyError:
            try:
                return self.index[obj._reversed.__uuid__] ^ 1
            except KeyError:
                raise KeyError(obj)

    def cache_all_cvs(self):
        for store, idx in self.cv_list.values():
            store.fill_cache()

    def pos(self, obj):
        return self.index.get(obj.__uuid__)

    def all(self):
        return peng.Trajectory(map(self.proxy, self.vars['uuid'][:]))

    def __getitem__(self, item):
        """
        Enable numpy style selection of object in the store
        """
        if type(item) is int or type(item) is str or type(item) is long:
            return self.load(item)
        elif type(item) is slice:
            return [self.load(idx)
                    for idx in range(*item.indices(len(self)))]
        elif type(item) is list:
            return [self.load(idx) for idx in item]
        elif item is Ellipsis:
            return iter(self)


class SnapshotValueStore(ObjectStore):
    def __init__(
            self,
            time_reversible=True,
            allow_incomplete=False,
            chunksize=250
    ):
        super(SnapshotValueStore, self).__init__(None)
        self.snapshot_index = None
        if not time_reversible and not allow_incomplete:
            raise RuntimeError(
                'Only time_reversible CVs can currently be '
                'stored using mode "complete"')

        self.time_reversible = time_reversible
        self.allow_incomplete = allow_incomplete
        self.chunksize = chunksize

        self.snapshot_pos = None
        self._len = 0

    def to_dict(self):
        return {
            'time_reversible': self.time_reversible,
            'allow_incomplete': self.allow_incomplete
        }

    def create_uuid_index(self):
        return dict()

    def register(self, storage, prefix):
        super(SnapshotValueStore, self).register(storage, prefix)
        self.snapshot_pos = self.storage.snapshots.pos

    def __len__(self):
        return len(self.variables['value'])

    # ==========================================================================
    # LOAD/SAVE DECORATORS FOR CACHE HANDLING
    # ==========================================================================

    def load(self, idx):
        pos = self.snapshot_pos(idx)

        if pos is None:
            return None

        if self.time_reversible:
            pos /= 2

        if self.allow_incomplete:
            # we want to load by uuid and it was not in cache.
            if pos in self.index:
                n_idx = self.index[pos]
            else:
                return None
            if n_idx < 0:
                return None
        else:
            if pos < self._len:
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

        if self.allow_incomplete:
            if pos in self.index:
                return

            n_idx = self.free()
            self.cache.update_size(n_idx)

        else:
            if pos < self._len:
                return

            n_idx = idx

        if self.allow_incomplete:
            # only if partial storage is used store index and update
            self.vars['index'][n_idx] = pos
            self.index[pos] = n_idx

        self.vars['value'][n_idx] = value
        self.cache[n_idx] = value
        self._len = max(self._len, n_idx + 1)

    def fill_cache(self):
        self.cache.load_max()

    def restore(self):
        if self.allow_incomplete:  # only if partial storage is used
            for pos, idx in enumerate(self.vars['index'][:]):
                self.index[idx] = pos

        self._len = len(self)
        self.initialize_cache()

    def initialize(self):
        self.initialize_cache()

    def initialize_cache(self):
        self.cache = LRUChunkLoadingCache(
            chunksize=self.chunksize,
            max_chunks=1000,
            variable=self.vars['value']
        )
        self.cache.update_size()

    def __getitem__(self, item):
        # enable numpy style selection of objects in the store
        try:
            if isinstance(item, peng.BaseSnapshot):
                return self.load(item)
            elif type(item) is list:
                return [self.load(idx) for idx in item]
            elif item is Ellipsis:
                return iter(self)
        except KeyError:
            return None

    def get(self, item):
        if self.allow_incomplete:
            try:
                return self[item]
            except KeyError:
                return None
        else:
            return self[item]
