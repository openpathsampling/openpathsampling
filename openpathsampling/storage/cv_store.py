from openpathsampling.netcdfplus import UniqueNamedObjectStore, NetCDFPlus, ObjectStore, VariableStore
from openpathsampling.netcdfplus.objects import UUIDDict

from openpathsampling.engines import BaseSnapshot

from uuid import UUID

class ObjectDictStore(UniqueNamedObjectStore):
    """
    ObjectStore to store a dict with StorableObject : value
    """
    def __init__(self, content_class, key_class):
        super(ObjectDictStore, self).__init__(
            content_class
        )
        self.key_class = key_class
        self._key_store = None

        self._cache_stores = dict()
        self._uuid_ref = None

    def initialize(self, units=None):
        super(ObjectDictStore, self).initialize()

        # index associated storage in class variable for all Trajectory instances to access

        if self.reference_by_uuid:
            self.storage.create_dimension('cv_cache')

            self.storage.create_variable(
                self.prefix + '_index',
                var_type='uuid',
                dimensions=('cv_cache',),
                description="the uuid for a specific cv_index",
                chunksizes=(10240,)
            )

    def to_dict(self):
        return {
            'content_class': self.content_class,
            'key_class': self.key_class
        }

    def register(self, storage, prefix):
        super(ObjectDictStore, self).register(storage, prefix)
        if self.reference_by_uuid:
            self._uuid_ref = UUIDDict()

    @property
    def key_store(self):
        """
        Return the associated store that contains the key elements

        Returns
        -------
        :class:`openpathsampling.netcdfplus.objects.ObjectStore`
            the object store instance
        """
        if self._key_store is None:
            self._key_store = self.storage.find_store(self.key_class)

        return self._key_store

    def _save(self, objectdict, idx):
        """
        Save the current state of the cache to the storage.

        Parameters
        ----------
        objectdict : :class:`openpathsampling.CollectiveVariable`
            the objectdict to store
        idx : int
            the index
        """
        self.vars['json'][idx] = objectdict

        if objectdict.cv_store_cache:
            self.create_cache(objectdict)

    def cache_var_name(self, idx):
        """
        Return the variable name use to store the values of the dict

        Parameters
        ----------
        idx : int
            the index of the objectdict of which the name is to be generated

        Returns
        -------
        str
            the name of the variable
        """
        if type(idx) is not int:
            idx = self.index.get(idx, None)
        if idx is not None:
            return 'cv%d' % idx

        raise KeyError("'%s' is neither an stored cv nor an integer index" % idx)

    def create_cache(self, cv):
        params = NetCDFPlus.get_value_parameters(cv(self.storage.template))
        shape = params['dimensions']

        if shape is None:
            chunksizes = None
        else:
            chunksizes = tuple(params['dimensions'])

        cache = KeyStore(cv)

        idx = self.index.get(cv, None)
        var_name = self.cache_var_name(idx)

        self.storage.register_store(var_name, cache)
        self.storage.create_dimension(cache.prefix, 0)

        if shape is not None:
            shape = tuple([0] + list(shape))
            chunksizes = tuple([1] + list(chunksizes))
        else:
            shape = tuple([0])
            chunksizes = tuple([1])

        cache.create_variable(
            'value',
            var_type=params['var_type'],
            dimensions=shape,
            chunksizes=chunksizes,
            simtk_unit=params['simtk_unit'],
        )

        cache.create_variable('index', 'index')

        print var_name

        self.set_cache_store(cv)

    def cache_var(self, obj):
        """
        Return the storage.vars[''] variable that contains the values

        Parameters
        ----------
        obj : :class:`openpathsampling.CollectiveVariable`
            the objectdict you request the attached variable store

        Returns
        -------
        :class:`openpathsampling.netcdfplus.netcdfplus.NetCDFPlus.ValueDelegate`

        """
        var_name = self.cache_var_name(obj)
        if var_name is None:
            return None

        return self.key_store.vars[var_name]

    def cache_store(self, cv):
        """
        Return the storage.vars[''] variable that contains the values

        Parameters
        ----------
        obj : :class:`openpathsampling.CollectiveVariable`
            the objectdict you request the attached variable store

        Returns
        -------
        :class:`openpathsampling.netcdfplus.ObjectStore` or `netcdf4.Variable`

        """
        idx = self.index.get(cv, None)
        var_name = self.cache_var_name(idx)

        return self.storage.stores[var_name]

    def has_cache(self, idx):
        return self.cache_var(idx) is not None

    def set_cache_store(self, objectdict):
        """
        Sets the attached storage of a CV to this store

        If you are using a CV in multiple files this allows to select which of the files is to be
        used as the file store. For performance reasons a single CV can be stored in multiple files, but its
        file cache can only be associated with a single store.

        This infers that if you have a full stored cache in one file and then save the CV in another file
        the old cache is still being used! If you switch to the new file you will have to recompute all CV
        values a second time. You can copy the store cache to the new file using `transfer_cache` but this
        requires that you have all snapshots loaded into memory otherwise there is no connection between
        snapshots. Meaning we will not try to figure out which pairs of snapshots in the two files are the
        same. You might have a snapshot stored in two files then remove the snapshot from memory and reload
        it. In this case the loaded snapshot will not know that it is also saved in another file.

        In the worst case you will have to compute the CVs again.

        Parameters
        ----------
        objectdict : :class:`openpathsampling.CollectiveVariable`
            the objectdict you want to set this store as its cache
        """
        idx = self.index.get(objectdict, None)
        if idx is not None:
            objectdict.set_cache_store(self.key_store, self.cache_var(idx))
        else:
            raise RuntimeWarning(('Your object is not stored as a CV in "%s" yet and hence a store ' +
                                  'for the cache cannot be attached.' +
                                 'Save your CV first and retry.') % self.storage)

    def cache_transfer(self, objectdict, target_file):
        """
        Transfer content of a stored objectdict from one file to another

        Parameters
        ----------
        objectdict : :class:`openpathsampling.CollectiveVariable`
            the objectdict you want to transfer
        target_file : objectdict : :class:`openpathsampling.netcdfplus.netcdfplus.NetCDFPlus`
            the target storage the cv should be transferred to

        """
        if objectdict in target_file.cvs.index:
            source_variable = self.cache_variable(objectdict)
            target_variable = target_file.cvs.cache_variable(objectdict)

            for source_idx, snapshot in enumerate(self.storage.snapshots):
                target_idx = target_file.snapshots.index.get(snapshot, None)
                if target_idx is not None:
                    target_variable[target_idx] = source_variable[source_idx]

    def sync(self, objectdict=None):
        """
        This will update the stored cache of the collective variable. It is
        different from saving in that the object is only created if it is
        saved (and the object caching will prevent additional creation)

        Parameters
        ----------
        objectdict : :class:`openpathsampling.CollectiveVariable` or `None`
            the objectdict to store. if `None` is given (default) then
            all collective variables are synced

        See also
        --------
        CollectiveVariable.sync

        """
        if objectdict is None:
            for cv in self:
                self.sync(cv)

            return

        objectdict.sync()

    def cache_all(self):
        """
        Fill the cache of all cvs

        """
        for cv in self:
            cv.cache_all()

    def _load(self, idx):
        # op = self.load_json(self.prefix + '_json', idx)
        op = self.vars['json'][idx]
        op.set_cache_store(self.key_store, self.cache_var(idx))

        return op

    def restore(self):
        if self.reference_by_uuid:
            self._uuid_ref = UUIDDict()

            for pos, idx in enumerate(self.vars['index'][:]):
                self._uuid_ref[self._uuid_ref[idx]] = pos

            self.load_indices()


class ReversibleObjectDictStore(ObjectDictStore):
    def cache_store(self, idx):

        var = self.cache_var(idx)
        if var is None:
            return None

        if var is not self._cache_stores:
            self._cache_stores[var] = self.storage.Key_Delegate(var, self.key_store)

        return self._cache_stores[var]

    def cache_transfer(self, objectdict, target_file):
        if objectdict in target_file.cvs.index:
            source_variable = self.cache_variable(objectdict)
            target_variable = target_file.cvs.cache_variable(objectdict)

            snapshots = objectdict.storage.snapshots.index.values()
            for snapshot in snapshots:
                source_idx = objectdict.storage.snapshots.index[snapshot]
                target_idx = target_file.snapshots.index.get(snapshot, None)
                if target_idx is not None:
                    target_variable[target_idx] = source_variable[source_idx]

    def cache_var(self, obj):
        var_name = self.cache_var_name(obj)
        if var_name is None:
            return None

        return self.storage._stores[var_name]

    def set_cache_store(self, objectdict):
        # print 'idx', self.index
        # print 'name', self.name_idx
        # print objectdict.__uuid__
        # print objectdict in self.index
        # print objectdict.__uuid__ in self.index
        # print self.index.get(objectdict)
        # print self.index.get(objectdict.__uuid__)
        # print self.idx(objectdict)

        idx = self.index.get(objectdict)
        if idx is not None:
            objectdict.set_cache_store(self.key_store, self.cache_var(idx))
        else:
            raise RuntimeWarning(('Your object is not stored as a CV in "%s" yet and hence a store ' +
                                  'for the cache cannot be attached.' +
                                 'Save your CV first and retry.') % self.storage)

    def _load(self, idx):
        op = self.vars['json'][idx]
        if op.cv_time_reversible:
            op.set_cache_store(self.key_store, self.cache_var(idx), self.cache_var(idx))
        else:
            op.set_cache_store(self.key_store, self.cache_var(idx), self.cache_bw(idx))

        return op

    def add_uuid(self, idx):
        if self.reference_by_uuid:
            length = int(len(self.variables['index']))
            print length, idx.__uuid__
            self._uuid_ref[idx] = length
            self.vars['index'][length] = idx.__uuid__

        return length


class KeyStore(ObjectStore):

    def __init__(self, cv):
        super(KeyStore, self).__init__(
            BaseSnapshot
        )

        self.cv = cv

    @property
    def uuid_ref(self):
        return self.storage.cvs._uuid_ref

    def create_int_index(self):
        return dict()

    # =============================================================================
    # LOAD/SAVE DECORATORS FOR CACHE HANDLING
    # =============================================================================

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

        # we want to load by uuid and it was not in cache.
        if idx in self.index:
            n_idx = self.index[idx]
        else:
            if self.fallback_store is not None:
                return self.fallback_store.load(idx)
            elif self.storage.fallback is not None:
                return self.storage.fallback.stores[self.name].load(idx)
            else:
                raise ValueError('str %s not found in storage or fallback' % idx)

        if n_idx < 0:
            return None

        # if it is in the cache, return it
        try:
            obj = self.cache[n_idx]
            return obj

        except KeyError:
            pass

        print n_idx

        obj = self.vars['value'][n_idx]

        self.index[idx] = n_idx
        self.cache[n_idx] = obj

        return obj

    def __setitem__(self, idx, value):
        """
        Saves an object to the storage.

        Parameters
        ----------
        obj : :py:class:`openpathsampling.netcdfplus.base.StorableObject`
            the object to be stored
        idx : int or string or `None`
            the index to be used for storing. This is highly discouraged since
            it changes an immutable object (at least in the storage). It is
            better to store also the new object and just ignore the
            previously stored one.

        """

        if idx in self.index:
            # has been saved so quit and do nothing
            return

        n_idx = self.free()
        self.vars['value'][n_idx] = value

        if idx in self.uuid_ref:
            self.vars['index'][n_idx] = self.uuid_ref[idx]
        else:
            self.vars['index'][n_idx] = self.storage.cvs.add_uuid(idx)

        self.index[idx] = n_idx
        self.cache[n_idx] = value

    def restore(self):
        uuid_ref = self.uuid_ref
        for pos, idx in enumerate(self.vars['index'][:]):
            self.index[uuid_ref[idx]] = pos

    def initialize(self):
        pass

    def __getitem__(self, item):
        """
        Enable numpy style selection of object in the store
        """
        try:
            if isinstance(item, BaseSnapshot):
                return self.load(item)
            elif type(item) is list:
                return [self.load(idx) for idx in item]
            elif item is Ellipsis:
                return self.iterator()
        except KeyError:
            return None
