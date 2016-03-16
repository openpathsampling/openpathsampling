from openpathsampling.netcdfplus import UniqueNamedObjectStore, NetCDFPlus


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

    def to_dict(self):
        return {
            'content_class': self.content_class,
            'key_class': self.key_class
        }

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
            return 'cv_%d_values' % idx

        raise KeyError("'%s' is neither an stored cv nor an integer index" % idx)

    def create_cache(self, objectdict):
        """
        Create the storage variable that holds the data for the object dict

        Parameters
        ----------
        objectdict : :class:`openpathsampling.CollectiveVariable`
            the object dictionary that you want the cache to be created for

        """
        idx = self.index.get(objectdict, None)
        if idx is not None:
            var_name = self.cache_var_name(idx)

            if var_name not in self.storage.variables:
                params = NetCDFPlus.get_value_parameters(objectdict(self.storage.template))

                shape = params['dimensions']

                if shape is None:
                    chunksizes = None
                else:
                    chunksizes = tuple(params['dimensions'])

                self.key_store.create_variable(
                    var_name,
                    var_type=params['var_type'],
                    dimensions=shape,
                    chunksizes=chunksizes,
                    simtk_unit=params['simtk_unit'],
                    maskable=True
                )
                self.storage.update_delegates()

        self.set_cache_store(objectdict)

    def cache_var(self, obj):
        """
        Return the storage.vars[''] variable that vontains the values

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

    def cache_variable(self, obj):
        """
        Return the storage.vars[''] variable that vontains the values

        Parameters
        ----------
        obj : :class:`openpathsampling.CollectiveVariable`
            the objectdict you request the attached variable store

        Returns
        -------
        :class:`openpathsampling.netcdf.ObjectStore` or `netcdf4.Variable`

        """
        var_name = self.cache_var_name(obj)
        snap_name = self.storage.snapshots.prefix

        return self.variables[snap_name + '_' + var_name]

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
        This will update the stored cache of the collectivevariable. It is
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


class ReversibleObjectDictStore(ObjectDictStore):

    def create_cache(self, objectdict):
        idx = self.index.get(objectdict, None)
        if idx is not None:
            var_name = self.cache_var_name(idx)

            if var_name not in self.storage.variables:

                params = NetCDFPlus.get_value_parameters(objectdict(self.storage.template))
                shape = params['dimensions']

                if shape is None:
                    chunksizes = None
                else:
                    chunksizes = tuple(shape)

                self.key_store.create_variable(
                    var_name + '_fw',
                    var_type=params['var_type'],
                    dimensions=shape,
                    chunksizes=chunksizes,
                    simtk_unit=params['simtk_unit'],
                    maskable=True
                )

                if not objectdict.cv_time_reversible:
                    self.key_store.create_variable(
                        var_name + '_bw',
                        var_type=params['var_type'],
                        dimensions=shape,
                        chunksizes=chunksizes,
                        simtk_unit=params['simtk_unit'],
                        maskable=True
                    )

                self.storage.update_delegates()

        self.set_cache_store(objectdict)

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

        return self.key_store.vars[var_name + '_fw']

    def cache_bw(self, obj):
        var_name = self.cache_var_name(obj)
        if var_name is None:
            return None

        return self.key_store.vars[var_name + '_bw']

    def set_cache_store(self, objectdict):
        idx = self.index.get(objectdict, None)
        if idx is not None:
            if objectdict.cv_time_reversible:
                objectdict.set_cache_store(self.key_store, self.cache_var(idx), self.cache_var(idx))
            else:
                objectdict.set_cache_store(self.key_store, self.cache_var(idx), self.cache_bw(idx))
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
