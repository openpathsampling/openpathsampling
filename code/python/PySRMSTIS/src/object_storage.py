import copy
import json
import yaml

def addcache(func):
    def inner(self, idx, *args, **kwargs):
        if idx in self.cache:
            return self.cache[idx]

        obj = func(self, idx, *args, **kwargs)
        self.cache[idx] = obj
        return obj
    return inner

class StoredObject(object):
    pass

class ObjectStorage(object):
    """
    Base Class for storing complex objects in a netCDF4 file. It holds a reference to the store file.
    """

    def __init__(self, storage, obj, named=False, json=False):
        """

        :param storage: a reference to the netCDF4 file
        :param obj: a reference to the Class to be stored
        :return:
        """
        self.storage = storage
        self.content_class = obj
        self.idx_dimension = obj.__name__.lower()
        self.cache = dict()
        self.named = named
        self.json = json

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
                self.cache[idx] = obj
                return None
            else:
                idx = self.free()
                obj.idx[storage] = idx
                self.cache[idx] = obj
                return idx
        else:
            idx = int(idx)
            self.cache[idx] = obj
            return idx

    def from_cache(self, idx):
        if idx in self.cache:
            return self.cache[idx]
        else:
            return None

    def to_cache(self, obj):
        self.cache[obj.idx[self.storage]] = obj

    @addcache
    def load(self, idx):
        '''
        Returns an object from the storage. Needs to be implented from the specific storage class.
        '''

        # Create object first to break any unwanted recursion in loading
        obj = StoredObject()
        setattr(obj, 'idx', {self.storage : idx})
        return self.load_json(self.idx_dimension + '_json', idx, obj)

    def save(self, obj, idx=None):
        '''
        Returns an object from the storage. Needs to be implented from the specific storage class.
        '''

        idx = self.index(obj, idx)

        if idx is not None:
            self.save_json(self.idx_dimension + '_json', idx, obj)


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
        print self.idx_dimension
        # define dimensions used for the specific object
        self.storage.createDimension(self.idx_dimension, 0)
        if self.named:
            self.init_variable(self.idx_dimension + "_name", 'str', description='A short descriptive name for convenience')
        if self.json:
            self.init_variable(self.idx_dimension + "_json", 'str', description='A json serialized version of the object')


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

    def save_mixed(self, name, idx, values, dimension=None):
        storage = self.storage
        begin = self.free_idx(name)
        length = len(values)
        for position, value in enumerate(values):
            storage.variables[name][begin + position] = value

        if dimension is None:
            dimension = name

        storage.variables[dimension + '_length'][idx] = length
        storage.variables[dimension + '_idx'][idx] = begin

    def load_mixed(self, name, idx, dimension=None):
        if dimension is None:
            dimension = name

        return self.storage.variables[name][self.slice(dimension, idx)]

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

    def save_object_list(self, name, idx, data, obj):
        specs = obj.details_specs
        self.storage.variables[name + '_reference_idx'][idx] = obj.idx[self.storage]
        self.save_list(name, idx, data, specs)

    def load_object_list(self, name, idx, cls):
        spec_obj_idx = self.storage.variables[name + '_reference_idx'][idx]
        store = getattr(self.storage, cls)
        spec_obj = store.load(spec_obj_idx)
        specs = spec_obj.specifications
        return self.load_list(name, idx, specs)

    def save_list(self, name, idx, data, specs):
        values = [None] * len(specs)

        for position, spec in enumerate(specs):
            if spec.cls is int:
                values[position] = float(data[spec.name])
            elif spec.cls is float:
                values[position] = float(data[spec.name])
            elif spec.cls is bool:
                values[position] = float(data[spec.name])
            elif spec.cls is str:
                # NOT IMPLEMENTED...
                pass
            else:
                index = data[spec.name].idx[self.storage]
                values[position] = float(index)

        self.save_mixed(name, idx, values)

    def load_list(self, name, idx, specs):
        ncfile = self.storage

        data = self.load_mixed(name, idx)

        if len(data) != len(specs):
            raise ValueError('Loaded data has different length of specification!')

        obj = {}
        for position, spec in enumerate(self.object_list_spec[name]):
            if spec.cls is int:
                r = int(data[position])
            elif spec.cls is float:
                r = float(data[position])
            elif spec.cls is bool:
                r = bool(data[position])
            elif spec.cls is str:
                # NOT IMPLEMENTED...
                pass
            else:
                store = getattr(ncfile, spec.cls)
                r = store.load(int(data[position]))

            obj[spec.name] = r

        return obj

    def init_list(self, name):
        ncfile = self.storage

        self.init_mixed_length(name)
        ncvar_name_idx  = ncfile.createVariable(name, 'f', name)

    def init_str(self, name):
        ncfile = self.storage

        self.init_mixed_length(name)
        ncvar_name_idx  = ncfile.createVariable(name, 'str', name)

    def _simplify_var(self,obj):
        if type(obj).__module__ != '__builtin__':
            if hasattr(obj, 'cls'):
                getattr(self.storage, obj.cls).save(obj)
                return { 'idx' : obj.idx[self.storage], 'cls' : obj.cls}
            else:
                return None
        elif type(obj) is list:
            return [self._simplify_var(o) for o in obj]
        elif type(obj) is dict:
            return {key : self._simplify_var(o) for key, o in obj.iteritems() if type(key) is str and key != 'idx' }
        else:
            return obj

    def _build_var(self,obj):
        if type(obj) is dict:
            if 'cls' in obj:
                return getattr(self.storage, obj['cls']).load(obj['idx'])
            else:
                return {key : self._build_var(o) for key, o in obj.iteritems()}
        elif type(obj) is list:
            return [self._build_var(o) for o in obj]
        else:
            return obj

    def save_json(self, name, idx, obj):
        data = obj.__dict__
        cls = obj.__class__.__name__
        store = obj.cls

        simplified = self._simplify_var([cls, store, data])
        json_string = json.dumps(simplified)

        self.storage.variables[name][idx] = json_string

    def load_json(self, name, idx, obj = None):
        json_string = self.storage.variables[name][idx]
        simplified = yaml.load(json_string)
        if obj is None:
            obj = StoredObject()

        data = self._build_var(simplified[2])
        for key, value in data.iteritems():
            setattr(obj, key, value)

        setattr(obj, 'cls', simplified[1])

        return obj


