from object_storage import ObjectStore

class ObjectVariableStore(ObjectStore):
    def __init__(self, cls, key_class, var_name):
        super(ObjectVariableStore, self).__init__(
            cls,
            has_uid=False,
            json=False,
            has_name=False
        )
        self.key_class = key_class
        self._key_store = None
        self.var_name = var_name

    @property
    def key_store(self):
        if self._key_store is None:
            self._key_store = self.storage._obj_store[self.key_class]

        return self._key_store

    def load(self, idx):
        return self.key_store.vars[idx]

    def save(self, idx, value):
        self.key_store.vars[idx] = value

class ObjectDictStore(ObjectStore):
    def __init__(self, cls, key_class):
        super(ObjectDictStore, self).__init__(
            cls,
            has_uid=True,
            json=True,
            has_name=True
        )
        self.key_class = key_class
        self._key_store = None

    @property
    def key_store(self):
        if self._key_store is None:
            self._key_store = self.storage._obj_store[self.key_class]

        return self._key_store

    def save(self, objectdict, idx=None):
        """
        Save the current state of the cache to the storage.

        Parameters
        ----------
        objectdict : object
            the objectdict to store
        idx : int
            the index
        """
        storage = self.storage

        if objectdict.store_cache:
            var_name = 'cv' + '_' + str(idx) + '_' + objectdict.name

            if var_name not in storage.variables:
                self.key_store.init_variable(
                    var_name,
                    objectdict.var_type,
                    objectdict.dimensions
                )
                self.storage.update_delegates()

        self.save_json(self.prefix + '_json', idx, objectdict)

        # this will copy the cache from an op and store it if it is stored
        if objectdict.store_cache:
            objectdict.flush_cache(storage)
            self.sync(objectdict)

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

        if objectdict.store_cache:
            objectdict.sync(store=self)

    def cache_all(self):
        for cv in self:
            cv.cache_all(self)

    def set_value(self, objectdict, position, value):
        idx = self.idx(objectdict)

        if idx is not None and idx >=0:
            var_name = 'cv_' + str(idx) + '_' + objectdict.name
            self.key_store.vars[var_name][position] = value

    def set_list_value(self, objectdict, positions, values):
        idx = self.idx(objectdict)

        if idx is not None and idx >=0:
            var_name = 'cv_' + str(idx) + '_' + objectdict.name
            self.key_store.vars[var_name][positions] = values

    def get_value(self, objectdict, position):
        idx = self.idx(objectdict)

        if idx is not None and idx >=0:
            var_name = self.key_store.prefix + '_cv_' + str(idx) + '_' + objectdict.name
            val = self.storage.variables[var_name][position]

            if hasattr(val, 'mask'):
                return None
            else:
                return self.storage.vars[var_name].getter(val)

        return None

    def get_list_value(self, objectdict, positions):
        idx = self.idx(objectdict)

        if idx is not None and idx >=0:
            var_name = self.key_store.prefix + '_cv_' + str(idx) + '_' + objectdict.name
            values = self.storage.variables[var_name][positions]
            getter = self.storage.vars[var_name].getter

            return [None if hasattr(val, 'mask') else getter(val) for val in values ]

        return [None] * len(positions)

    def load(self, idx):
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

        op = self.load_json(self.prefix + '_json', idx)

        return op

    def _init(self, **kwargs):
        """
        Initialize the associated storage to allow for ensemble storage

        """
        super(ObjectDictStore, self)._init()