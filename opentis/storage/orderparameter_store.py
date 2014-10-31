from object_storage import ObjectStorage
from wrapper import loadcache, savecache

class ObjectDictStorage(ObjectStorage):

    def __init__(self, storage, cls, value_type):
        super(ObjectDictStorage, self).__init__(storage, cls, named=True)
        self.idx_dimension = 'dict_' + self.idx_dimension
        self.value_type = value_type

    def save(self, objectdict, idx=None):
        """
        Save the current state of the cache to the storage.

        Parameters
        ----------
        storage : Storage() on None
            The storage (not ObjectStorage) to store in. If None then all associated storages will be saved in.

        """
        storage = self.storage
        if idx is None:
            if storage in objectdict.idx:
                self.cache[idx] = objectdict
                idx = objectdict.idx[storage]
            else:
                idx = self.free()
                objectdict.idx[storage] = idx
                self.cache[idx] = objectdict
        else:
            idx = int(idx)
            self.cache[idx] = objectdict

        self._update_store(objectdict)
        store = objectdict.storage_caches[storage]
        self.save_objectdict(self.idx_dimension, int(idx), store, self.value_type)
        self.tidy_cache(objectdict)

    @loadcache
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
        data = self.load_objectdict(self.idx_dimension,int(idx), self.content_class.__name__.lower(), self.value_type)
        if op is None:
            # create StorageObject
            pass

        op.storage_caches[storage] = data

    def _init(self):
        """
        Initialize the associated storage to allow for ensemble storage

        """
        super(ObjectDictStorage, self)._init()

        self.init_objectdict(self.idx_dimension, self.content_class.__name__.lower())

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