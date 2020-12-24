import logging
from uuid import UUID

import openpathsampling.engines as peng
from openpathsampling.netcdfplus import ObjectStore, \
    NetCDFPlus, LoaderProxy

from .snapshot_feature import FeatureSnapshotStore
from .snapshot_value import SnapshotValueStore

import sys
if sys.version_info > (3, ):
    unicode = str
    long = int

logger = logging.getLogger(__name__)
init_log = logging.getLogger('openpathsampling.initialization')


class ReversalHashedList(dict):
    def __init__(self):
        dict.__init__(self)
        self._list = []

    def append(self, key):
        dict.__setitem__(self, key & ~1, len(self._list) * 2 ^ (key & 1))
        self._list.append(key)

    # noinspection PyCallByClass
    def extend(self, t):
        l = len(self._list)
        # t = filter(t, lambda x : x not in self)
        dict.update(
            self,
            map(lambda x, y: (x & ~1, y * 2 ^ (x & 1)), t, range(l, l + len(t)))
        )
        # map(lambda x, y: dict.__setitem__(self, x & ~1, y * 2 ^ (x & 1)), t,
        #     range(l, l + len(t)))
        self._list.extend(t)

    def __len__(self):
        return len(self._list) * 2

    def __setitem__(self, key, value):
        # we will always store the ones with even keys
        dict.__setitem__(self, key & ~1, value ^ (key & 1))
        # we will always store the ones with even value
        self._list[value // 2] = key ^ (value & 1)

    def get(self, key, d=None):
        uu = dict.get(self, key & ~1, d)
        if uu is not None:
            return uu ^ (key & 1)
        else:
            return d

    def __getitem__(self, key):
        return dict.__getitem__(self, key & ~1) ^ (key & 1)

    def __contains__(self, key):
        return dict.__contains__(self, key & ~1)

    def index(self, key):
        return self._list[key // 2] ^ (key & 1)

    def mark(self, key):
        k = key & ~1
        if k not in self:
            dict.__setitem__(self, k, -2)

    def unmark(self, key):
        k = key & ~1
        if k in self:
            dict.__delitem__(self, k)

    @property
    def list(self):
        return self._list


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
        self._store = {}

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

        if isinstance(idx, (int, long)):
            # we want to load by uuid and it was not in cache.
            if idx < 10000000000:
                n_idx = idx
            elif idx in self.index:
                n_idx = self.index[idx]
            else:
                if self.fallback_store is not None:
                    return self.fallback_store.load(idx)
                elif self.storage.fallback is not None:
                    return self.storage.fallback.stores[self.name].load(idx)
                else:
                    raise ValueError(
                        'int %s not found in storage or fallback' % idx)

        elif type(idx) is not int:
            raise ValueError(
                ('indices of type "%s" are not allowed in named storage '
                 '(only str and int)') % type(idx).__name__
            )
        else:
            n_idx = int(idx)

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
        store_idx = int(self.variables['store'][idx // 2])

        if store_idx < 0:
            if self.fallback_store is not None:
                return self.fallback_store.load(idx)
            elif self.storage.fallback is not None:
                return self.storage.fallback.snapshots.load(idx)
            else:
                raise KeyError(
                    '%s not found in storage or fallback' % idx)
        else:
            store = self.store_snapshot_list[store_idx]
            snap = store[int(idx)]
            return snap

    def __len__(self):
        return len(self.storage.dimensions[self.prefix]) * 2

    def initialize(self):
        super(SnapshotWrapperStore, self).initialize()

        self.create_variable('store', 'index')

        self.storage.create_dimension('snapshottype')

        self.storage.create_variable(
            'snapshottype',
            'obj.stores',
            'snapshottype')

    def free(self):
        idx = len(self)
        while idx in self._free:
            idx += 2

        return idx

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
        super(SnapshotWrapperStore, self).restore()

        for idx, store in enumerate(self.storage.vars['snapshottype']):
            self.type_list[store.descriptor] = (store, idx)
            self.store_snapshot_list.append(store)

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
        n_idx = self.index.get(obj.__uuid__)

        if n_idx is not None:
            # snapshot is mentioned
            store_idx = int(self.variables['store'][n_idx // 2])
            if not store_idx == -1:
                # and stored
                return self.reference(obj)

        if self.only_mention:
            if n_idx is None:
                n_idx = len(self.index)

                # only mention but not really store snapshots
                self.vars['store'][n_idx // 2] = -1
                self.index.append(obj.__uuid__)
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

        if isinstance(obj, LoaderProxy):
            if obj._store is self:
                # is a proxy of a saved object so do nothing
                return obj.__uuid__
            else:
                # it is stored but not in this store so we try storing the
                # full attribute which might be still in cache or memory
                # if that is not the case it will be stored again. This can
                # happen when you load from one store save to another. And load
                # again after some time while the cache has been changed and try
                # to save again the loaded object. We will not explicitly store
                # a table that matches objects between different storages.
                return self.save(obj.__subject__)

        if n_idx is None:
            n_idx = len(self.index)

            self._save(obj, n_idx)
            self._auto_complete_single_snapshot(obj, n_idx)
            self._set_id(n_idx, obj)
        else:
            store, store_idx = self.type_list[obj.engine.descriptor]
            self.vars['store'][n_idx // 2] = store_idx
            store[n_idx] = obj

        self.cache[n_idx] = obj

        return self.reference(obj)

    def _save(self, obj, n_idx):
        try:
            store, store_idx = self.type_list[obj.engine.descriptor]
            self.vars['store'][n_idx // 2] = store_idx
            self.index.append(obj.__uuid__)
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
                self.vars['store'][n_idx // 2] = store_idx
                self.index.append(obj.__uuid__)
                store[n_idx] = obj
                return store

            elif self.treat_missing_snapshot_type == 'ignore':
                # we keep silent about it
                self.vars['store'][n_idx // 2] = -1
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
        for cv, cv_store in self.attribute_list.items():
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
                            n_idx = pos // 2
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
        if cv not in self.attribute_list:
            return

        cv_store = self.attribute_list[cv]

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

        if cv not in self.attribute_list:
            return

        cv_store = self.attribute_list[cv]

        # for complete this does not make sense
        if cv_store.allow_incomplete:

            # loop all objects in the fast CV cache
            for obj, value in cv._cache_dict.cache.items():
                if value is not None:
                    pos = self.index.get(obj.__uuid__)

                    # if the snapshot is not saved, nothing we can do
                    if pos is None:
                        continue

                    if cv_store.time_reversible:
                        pos //= 2

                    if pos in cv_store.index:
                        # this value is stored so skip it
                        continue

                    n_idx = cv_store.free()

                    cv_store.vars['value'][n_idx] = value
                    cv_store.vars['index'][n_idx] = pos
                    cv_store.index[pos] = n_idx
                    cv_store.cache[n_idx] = value

    @staticmethod
    def _get_cv_name(cv_idx):
        return 'cv' + str(cv_idx)

    def add_attribute(
            self, store_cls, attribute, template,
            allow_incomplete=None, chunksize=None):

        self.add_cv(attribute, template, allow_incomplete, chunksize)

    # todo: this can be reduced to almost the function in super
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
        if cv in self.attribute_list:
            return self.attribute_list[cv]

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

            # todo: Remove once netcdf4/python bug in Unidata/netcdf4-python#566 has
            # been resolved. It seems that so far, after data has been written to a
            # dimension you have to use chunksize 1 for any new variables in that
            # dimension Even if the other variables all share the same chunksize
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

        store.initialize()

        self.attribute_list[cv] = store
        attribute_idx = self.storage.cvs.index[cv.__uuid__]
        self.storage.attributes.vars['cache'][attribute_idx] = store

        # use the cache and function of the CV to fill the store when it is made
        if not allow_incomplete:

            indices = self.vars['uuid'][:]

            for pos, idx in enumerate(indices):

                proxy = LoaderProxy.new(self.storage.snapshots, idx)
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
        return store

    def create_uuid_index(self):
        return ReversalHashedList()

    def _get_id(self, idx, obj):
        uuid = self.index.index(int(idx))
        obj.__uuid__ = uuid
        if obj._reversed:
            obj._reversed.__uuid__ = uuid ^ 1

    def _set_id(self, idx, obj):
        self.vars['uuid'][idx // 2] = obj.__uuid__

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

    def all(self):
        return peng.Trajectory(map(self.proxy, self.index.list))

    # # todo: this will not catch a non found item as the super function does!
    # def __getitem__(self, item):
    #     """
    #     Enable numpy style selection of object in the store
    #     """
    #     if type(item) is int:
    #         if item < 0:
    #             item += len(self)
    #         return self.load(item)
    #     elif type(item) is str or type(item) is long:
    #         return self.load(item)
    #     elif type(item) is slice:
    #         return [self.load(idx)
    #                 for idx in range(*item.indices(len(self)))]
    #     elif type(item) is list:
    #         return [self.load(idx) for idx in item]
    #     elif item is Ellipsis:
    #         return iter(self)
