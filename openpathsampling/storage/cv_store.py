from openpathsampling.netcdfplus import ObjectStore


class ObjectDictStore(ObjectStore):
    def __init__(self, content_class, key_class):
        super(ObjectDictStore, self).__init__(
            content_class,
            json=True,
            has_name=True
        )
        self.key_class = key_class
        self._key_store = None

        self._cache_stores = dict()

    def to_dict(self):
        return {
            'content_class': self.content_class,
            'key_class': self.key_class
        }

    @property
    def key_store(self):
        if self._key_store is None:
            self._key_store = self.storage._obj_store[self.key_class]

        return self._key_store

    def _save(self, objectdict, idx):
        """
        Save the current state of the cache to the storage.

        Parameters
        ----------
        objectdict : object
            the objectdict to store
        idx : int
            the index
        """
        self.vars['json'][idx] = objectdict

        if objectdict.store_cache:
            self.create_cache(objectdict)


    def cache_var_name(self, idx):
        if type(idx) is not int:
            idx = self.index.get(idx, None)
        if idx is not None:
            return 'cv_%d_values' % idx

        raise KeyError("'%s' is neither an stored cv nor an integer index" % idx)

    def create_cache(self, objectdict):
        idx = self.index.get(objectdict, None)
        if idx is not None:
            var_name = self.cache_var_name(idx)

            if var_name not in self.storage.variables:

                params = objectdict.return_parameters_from_template(self.storage.template)

                self.key_store.init_variable(
                    var_name,
                    var_type=params['cv_return_type'],
                    dimensions=params['cv_return_shape'],
                    simtk_unit=params['cv_return_simtk_unit'],
                    maskable=True
                )
                self.storage.update_delegates()

        self.set_cache_store(objectdict)

    def cache_var(self, obj):
        var_name = self.cache_var_name(obj)
        if var_name is None:
            return None

        return self.key_store.vars[var_name]

    def cache_variable(self, obj):
        var_name = self.cache_var_name(obj)
        snap_name = self.storage.snapshots.prefix

        return self.variables[snap_name + '_' + var_name]

    def cache_store(self, idx):

        var = self.cache_var(idx)
        if var is None:
            return None

        if var is not self._cache_stores:
            self._cache_stores[var] = self.storage.Key_Delegate(var, self.key_store)

        return self._cache_stores[var]

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
        """
        idx = self.index.get(objectdict, None)
        if idx is not None:
            objectdict.set_cache_store(self.key_store, self.cache_var(idx))
        else:
            raise RuntimeWarning(('Your object is not stored as a CV in "%s" yet and hence a store ' +
                                  'for the cache cannot be attached.' +
                                 'Save your CV first and retry.') % self.storage)

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

    def sync(self, objectdict=None):
        """
        This will update the stored cache of the collectivevariable. It is
        different from saving in that the object is only created if it is
        saved (and the object caching will prevent additional creation)

        Parameters
        ----------
        objectdict : object or None (default)
            the objectdict to store. if `None` is given (default) then
            all collectivevariables are synced

        See also
        --------
        CollectiveVariable.sync

        """
        if objectdict is None:
            for obj in self:
                self.sync(obj)

            return

        objectdict.sync()

    def cache_all(self):
        for cv in self:
            cv.cache_all()

    def _load(self, idx):
        """
        Restores the cache from the storage using the name of the
        collectivevariable.

        Parameters
        ----------
        storage : Storage() on None
            The storage (not ObjectStore) to store in. If None then all
            associated storages will be loaded from.

        Notes
        -----
        Make sure that you use unique names otherwise you might load the
        wrong parameters!
        """

        # op = self.load_json(self.prefix + '_json', idx)
        op = self.vars['json'][idx]
        op.set_cache_store(self.key_store, self.cache_var(idx))

        return op

    def _init(self, **kwargs):
        """
        Initialize the associated storage to allow for ensemble storage

        """
        super(ObjectDictStore, self)._init()
