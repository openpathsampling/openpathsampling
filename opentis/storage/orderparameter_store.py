from object_storage import ObjectStorage
from wrapper import loadcache, savecache

class ObjectDictStorage(ObjectStorage):

    def __init__(self, storage, cls, value_type):
        super(ObjectDictStorage, self).__init__(storage, cls, named=True)
        self.idx_dimension = 'dict_' + self.idx_dimension
        self.value_type = value_type

    @savecache
    def save(self, objectdict, idx=None):
        """
        Save the current state of the cache to the storage.

        Parameters
        ----------
        storage : Storage() on None
            The storage (not ObjectStorage) to store in. If None then all associated storages will be saved in.

        """
        storage = self.storage

        objectdict._update_store(storage)
        store = objectdict.storage_caches[storage]
        self.save_objectdict(self.idx_dimension, int(idx), store, self.value_type)


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

    def _update_store(self, storage = None):
        """
        This will transfer everything from the memory cache into the storage copy in memory which is used to interact with
        the file storage.

        Parameters
        ----------
        storage : Storage() on None
            The storage (not ObjectStorage) to store in. If None then all associated storages will be updated up.

        """
        if storage is None:
            if len(self.storage_caches) > 0:
                map(self._update_store, self.storage_caches.keys())
        else:
            if storage not in self.storage_caches:
                # TODO: Throw exception
                self.storage_caches[storage] = dict()

            store = self.storage_caches[storage]
            for item, value in self.iteritems():
                if storage in item.idx:
                    store[item.idx[storage]] = value

    def tidy_cache(self, storage = None):
        """
        This will transfer everything from the memory cache into the storage copy in memory which is used to interact with
        the file storage.

        Parameters
        ----------
        storage : Storage() on None
            The storage (not ObjectStorage) to store in. If None then all associated storages will be cleaned up.

        """

        # TODO: This doesnt work because the storage is changed during deleting superfluous elements
        if storage is None:
            if len(self.storage_caches) > 0:
                map(self.tidy_cache, self.storage_caches.keys())
        else:
            # Make sure configuration_indices are stored and have an index and then add the configuration index to the trajectory

            if storage not in self.storage_caches:
                # TODO: Throw exception
                self.storage_caches[storage] = dict()

            new_dict = {item:value for item, value in self.iteritems() if storage not in item.idx}

            self.clear()
            self.update(new_dict)