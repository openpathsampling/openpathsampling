from object_storage import ObjectStorage
from opentis.orderparameter import OrderParameter

class ObjectDictStorage(ObjectStorage):

    def __init__(self, storage, cls, key_class):
        super(ObjectDictStorage, self).__init__(storage, cls, named=True, identifier='name', load_lazy=False)
        self.key_class = key_class

    def save(self, objectdict, idx=None):
        """
        Save the current state of the cache to the storage.

        Parameters
        ----------
        """
        storage = self.storage

        self._update_store(objectdict)
        store = objectdict.storage_caches[storage]
        length = len(store)

        var_name = self.idx_dimension + '_' + str(idx) + '_' + objectdict.name

        if var_name + '_value' not in self.storage.variables:
            self.init_variable(var_name + '_value', 'float', (self.key_class.__name__.lower()))
            self.init_variable(var_name + '_set', 'index', (self.key_class.__name__.lower()))

        self.storage.variables[self.idx_dimension + '_name'][idx] = objectdict.name
        self.save_variable(self.idx_dimension + '_length', idx, length)
        self.storage.variables[var_name + '_value'][store.keys()] = self.list_to_numpy(store.values(), 'float')
        self.storage.variables[var_name + '_set'][0:length] = store.keys()

        self.tidy_cache(objectdict)

    def load(self, idx, op=None):
        """
        Restores the cache from the storage using the name of the orderparameter.

        Parameters
        ----------
        storage : Storage() on None
            The storage (not ObjectStorage) to store in. If None then all associated storages will be loaded from.

        Notes
        -----
        Make sure that you use unique names otherwise you might load the wrong parameters!
        """

        storage = self.storage
#        data = self.load_objectdict(self.idx_dimension,int(idx), self.content_class.__name__.lower(), float)

        name = storage.variables[self.idx_dimension + '_name'][idx]
        var_name = self.idx_dimension + '_' + str(idx) + '_' + name
        length = self.load_variable(self.idx_dimension + '_length', idx)
        stored_idx = self.storage.variables[var_name + '_set'][0:length]
        data_all = storage.variables[var_name + '_value'][:]
        data = self.list_from_numpy(data_all[self.list_from_numpy(stored_idx, 'index')], 'float')

        if op is None:
            op = OrderParameter(name)

        op.storage_caches[storage] = dict(zip(stored_idx, data))

        return op

    def restore(self, obj):
        idx = self.idx_by_name(obj.identifier)

        if idx is not None:
            return self.load(idx, obj)
        else:
            return None

    def _init(self):
        """
        Initialize the associated storage to allow for ensemble storage

        """
        super(ObjectDictStorage, self)._init()

        self.init_variable(self.idx_dimension + '_length', 'index', self.idx_dimension, chunksizes=(1, ))

    def _update_store(self, obj):
        """
        This will transfer everything from the memory cache into the storage copy in memory which is used to interact with
        the file storage.

        Parameters
        ----------
        storage : Storage() on None
            The storage (not ObjectStorage) to store in. If None then all associated storages will be updated up.

        """

        storage = self.storage

        if storage not in obj.storage_caches:
            # TODO: Throw exception
            obj.storage_caches[storage] = dict()

        store = obj.storage_caches[storage]
        for item, value in obj.iteritems():
            if storage in item.idx:
                store[item.idx[storage]] = value

    def tidy_cache(self, obj):
        """
        This will transfer everything from the memory cache into the storage copy in memory which is used to interact with
        the file storage.

        Parameters
        ----------
        storage : Storage() on None
            The storage (not ObjectStorage) to store in. If None then all associated storages will be cleaned up.

        """

        storage = self.storage

        if storage not in obj.storage_caches:
            # TODO: Throw exception
            obj.storage_caches[storage] = dict()

        new_dict = {item: value for item, value in obj.iteritems() if storage not in item.idx}

        obj.clear()
        obj.update(new_dict)
