from object_storage import ObjectStore
from openpathsampling.collectivevariable import CollectiveVariable

class ObjectDictStore(ObjectStore):
    def __init__(self, storage, cls, key_class):
        super(ObjectDictStore, self).__init__(storage, cls, has_uid=True, json=False)
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
        storage = self.storage

        var_name = self.idx_dimension + '_' + str(idx) + '_' + objectdict.name

        if var_name + '_value' not in self.storage.variables:
            self.init_variable(var_name + '_value', 'float', (self.key_class.__name__.lower()))

        storage.variables[self.idx_dimension + '_name'][idx] = objectdict.name

        # this will copy the cache from an op and store it
        objectdict.flush_cache(self.storage)
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

        objectdict.sync(storage=self.storage)
#        self.storage.sync()

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

            if hasattr(val, 'mask'):
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

        storage = self.storage

        name = storage.variables[self.idx_dimension + '_name'][idx]
        op = CollectiveVariable(name)

        return op

    def _init(self):
        """
        Initialize the associated storage to allow for ensemble storage

        """
        super(ObjectDictStore, self)._init()

        self.init_variable(
            self.idx_dimension + '_name',
            'str',
            self.idx_dimension, chunksizes=(1, )
        )
