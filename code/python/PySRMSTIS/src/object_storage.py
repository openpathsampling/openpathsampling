import copy

class ObjectStorage(object):
    """
    Base Class for storing complex objects in a netCDF4 file. It holds a reference to the store file.
    """

    def __init__(self, storage, obj, register=False):
        """

        :param storage: a reference to the netCDF4 file
        :param obj: a reference to the Class to be stored
        :return:
        """
        self.storage = storage
        self.content_class = obj
        self.idx_dimension = obj.__name__.lower()
        self.cache = dict()

    def register(self):
        self.storage.links.append(self)
        return self

    def copy(self):
        """
        Create a deep copy of the ObjectStorage instance

        Returns
        -------
        ObjectStorage()
            the copied instance
        """
        store = copy.deepcopy(self)
        return store

    def __call__(self, storage):
        """
        Create a deep copy of the ObjectStorage instance using the new store provided as function argument

        Returns
        -------
        ObjectStorage()
            the copied instance
        """
        store = self.copy()
        store.storage = storage
        return store

    def index(self, obj, idx=None):
        """
        Return the appropriate index for saving or None if the object is already stored!

        Returns
        -------
        int
            the index that should be used or None if no saving is required.
        """

        storage = self.storage
        if idx is None:
            if storage in obj.idx:
                # has been saved so quit and do nothing
                return None
            else:
                idx = self.free()
                obj.idx[storage] = idx
                return idx
        else:
            return idx

    def from_cache(self, idx):
        if idx in self.cache:
            return self.cache[idx]
        else:
            return None

    def to_cache(self, obj):
        self.cache[obj.idx[self.storage]] = obj

    def load(self, idx):
        '''
        Returns an object from the storage. Needs to be implented from the specific storage class.
        '''
        return None

    def get(self, indices):
        """
        Returns a list of objects from the given list of indices

        Arguments
        ---------
        indices : list of int
            the list of integers specifying the object to be returned

        Returns
        -------
        list of objects
            a list of objects stored under the given indices

        """
        return [self.load(idx) for idx in range(0, self.count())[indices] ]


    def last(self):
        '''
        Returns the last generated trajectory. Useful to continue a run.

        Returns
        -------
        Trajectoy
            the actual trajectory object
        '''
        return self.load(self.count())

    def first(self):
        '''
        Returns the last stored object. Useful to continue a run.

        Returns
        -------
        Object
            the actual last stored object
        '''
        return self.load(1)

    def count(self):
        '''
        Return the number of objects in the storage

        Returns
        -------
        number : int
            number of objects in the storage.
        '''
        length = int(len(self.storage.dimensions[self.idx_dimension])) - 1
        if length < 0:
            length = 0
        return length


    def free(self):
        '''
        Return the number of the next free index

        Returns
        -------
        index : int
            the number of the next free index in the storage. Used to store a new object.
        '''
        return  self.count() + 1

    def _init(self):
        """
        Initialize the associated storage to allow for object storage. Mainly creates an index dimension with the name of the object.

        """
        # define dimensions used for the specific object
        self.storage.createDimension(self.idx_dimension, 0)

    def slice(self, name, idx):
        begin = self.idx(name, idx)
        end = begin + self.len(name, idx)
        return slice(begin, end)

    def idx(self, name, idx):
        return int(self.storage.variables[name + '_idx'][int(idx)])

    def len(self, name, idx):
        return int(self.storage.variables[name + '_length'][int(idx)])

    def free_idx(self, name):
        length = int(len(self.storage.dimensions[name]))
        return length

    def save_mixed(self, name, idx, values):
        storage = self.storage
        begin = self.free_idx(name)
        length = len(values)
        for position, value in enumerate(values):
            storage.variables[name][begin + position] = value

        storage.variables[name + '_length'][idx] = length
        storage.variables[name + '_idx'][idx] = begin

    def init_mixed_length(self, name, dimension=None):
        # index associated storage in class variable for all Origin instances to access
        ncfile = self.storage

        # define dimensions used in trajectories
        ncfile.createDimension(name, 0)                     # unlimited number of iterations

        if dimension is None:
            dimension = self.idx_dimension

        # Create variables for trajectories
        ncvar_name_idx  = ncfile.createVariable(name + '_idx', 'u4', dimension)
        ncvar_name_length  = ncfile.createVariable(name + '_length', 'u4', dimension)

        # Define units for snapshot variables.
        setattr(ncvar_name_idx,      'units', 'none')
        setattr(ncvar_name_length,      'units', 'none')

    def init_variable(self, name, var_type, dimensions = None, units='none', description=None):

        ncfile = self.storage

        if dimensions is None:
            dimensions = self.idx_dimension

        ncvar     = ncfile.createVariable(name, var_type, dimensions)

        # Define units for snapshot variables.
        setattr(ncvar,      'units', units)

        if description is not None:
            # Define long (human-readable) names for variables.
            setattr(ncvar,    "long_str", description)

    def save_object(self, var, idx, obj):
        self.storage.variables[var + '_idx'][idx] = obj.idx[self.storage]
