from object_storage import ObjectStore
from openpathsampling.orderparameter import OrderParameter

import numpy as np

class ObjectDictStore(ObjectStore):
    def __init__(self, storage, cls, key_class):
        super(ObjectDictStore, self).__init__(storage, cls, is_named=True, json=False)
        self.key_class = key_class

    def save(self, objectdict, idx):
        """
        Save the current state of the cache to the storage.

        Parameters
        ----------
        objectdict : object
            the objectdict to store
        idx : int
            the index
        """
        var_name = self.idx_dimension + '_' + str(idx) + '_' + objectdict.name

        if var_name + '_value' not in self.storage.variables:
            self.init_variable(var_name + '_value', 'float', (self.key_class.__name__.lower()))

        self.sync(objectdict)

    def sync(self, objectdict=None):
        """
        This will update the stored cache of the orderparameter. It is
        different from saving in that the object is only created if it is
        saved (and the object caching will prevent additional creation)

        Parameters
        ----------
        objectdict : object or None (default)
            the objectdict to store. if `None` is given (default) then
            all orderparameters are synced

        """
        storage = self.storage
        idx = self.idx(objectdict)

        if idx is not None and idx >=0:
            cache = objectdict.store_dict
#            length = len(cache)

            var_name = self.idx_dimension + '_' + str(idx) + '_' + objectdict.name

#            self.save_variable(self.idx_dimension + '_length', idx, length)
            storage.variables[var_name + '_value'][cache.keys()] = self.list_to_numpy(cache.values(), 'float')

    def set_value(self, objectdict, position, value):
        storage = self.storage
        idx = self.idx(objectdict)

        if idx is not None and idx >=0:
            var_name = self.idx_dimension + '_' + str(idx) + '_' + objectdict.name
            storage.variables[var_name + '_value'][position] = value

    def set_list_value(self, objectdict, positions, values):
        storage = self.storage
        idx = self.idx(objectdict)

        if idx is not None and idx >=0:
            var_name = self.idx_dimension + '_' + str(idx) + '_' + objectdict.name
            storage.variables[var_name + '_value'][positions] = values


    def get_value(self, objectdict, position):
        storage = self.storage
        idx = self.idx(objectdict)

        if idx is not None and idx >=0:
            var_name = self.idx_dimension + '_' + str(idx) + '_' + objectdict.name
            val = storage.variables[var_name + '_value'][position]

            if val is np.ma.masked:
                return None
            else:
                return val

        return None

    def get_list_value(self, objectdict, positions):
        storage = self.storage
        idx = self.idx(objectdict)

        if idx is not None and idx >=0:
            var_name = self.idx_dimension + '_' + str(idx) + '_' + objectdict.name
            val = storage.variables[var_name + '_value'][positions]

            return val.tolist()

        return [None] * len(positions)


    def load(self, idx):
        """
        Restores the cache from the storage using the name of the
        orderparameter.

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

        storage = self.storage

        name = storage.variables[self.idx_dimension + '_name'][idx]
        var_name = self.idx_dimension + '_' + str(idx) + '_' + name

#        length = self.load_variable(self.idx_dimension + '_length', idx)
#        stored_idx = storage.variables[var_name + '_set'][0:length]
#        data_all = storage.variables[var_name + '_value'][:]
#        data = self.list_from_numpy(data_all[self.list_from_numpy(stored_idx, 'index')], 'float')

        op = OrderParameter(name)

        return op

    def _init(self):
        """
        Initialize the associated storage to allow for ensemble storage

        """
        super(ObjectDictStore, self)._init()

#        self.init_variable(self.idx_dimension + '_length', 'index', self.idx_dimension, chunksizes=(1, ))

    # def _update_store(self, obj):
    #     """
    #     This will transfer everything from the memory cache into the storage
    #     copy in memory which is used to interact with the file storage.
    #
    #     Parameters
    #     ----------
    #     storage : Storage() on None
    #         The storage (not ObjectStore) to store in. If None then all
    #         associated storages will be updated up.
    #
    #     """
    #
    #     storage = self.storage
    #
    #     if storage not in obj.storage_caches:
    #         # TODO: Throw exception
    #         obj.storage_caches[storage] = dict()
    #
    #     store = obj.storage_caches[storage]
    #     for item, value in obj.iteritems():
    #         if storage in item.idx:
    #             store[item.idx[storage]] = value
    #
    # def tidy_cache(self, obj):
    #     """
    #     This will transfer everything from the memory cache into the storage copy in memory which is used to interact with
    #     the file storage.
    #
    #     Parameters
    #     ----------
    #     storage : Storage() on None
    #         The storage (not ObjectStore) to store in. If None then all associated storages will be cleaned up.
    #
    #     """
    #
    #     storage = self.storage
    #
    #     if storage not in obj.storage_caches:
    #         # TODO: Throw exception
    #         obj.storage_caches[storage] = dict()
    #
    #     new_dict = {item: value for item, value in obj.iteritems() if storage not in item.idx}
    #
    #     obj.clear()
    #     obj.update(new_dict)