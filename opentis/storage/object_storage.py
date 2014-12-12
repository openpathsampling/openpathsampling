import copy

import yaml
import numpy as np

from wrapper import savecache, saveidentifiable, loadcache
from util import StorableObjectJSON
import simtk.unit as u

def add_storage_name(func):
    def inner(self, name, *args, **kwargs):
        var_name = '_'.join([self.db, name])
        func(self, var_name, *args, **kwargs)

    return inner

class StoredObject(object):
    def __getattr__(self, item):
        if hasattr(self, '_loader'):
            self._loader()
            delattr(self, '_loader')

        return self.__dict__[item]

    pass

# TODO: Combine the cache and all_names to be stored in one bis dict

class ObjectStorage(object):
    """
    Base Class for storing complex objects in a netCDF4 file. It holds a reference to the store file.
    """

    def __init__(self, storage, obj, named=False, json=False, identifier=None, dimension_units=None):
        """
        Attributes
        ----------

        storage : Storage
        content_class : class
            a reference to the class type to be stored using this Storage
        idx_dimension : str
            name of the dimension used for major numbering the stored object in the netCDF file.
            This is usually the lowercase class name
        cache : dict {int : object}
            A dictionary pointing to the actual stored object by index. It is filled at saving or loading
        named : bool
            Set, if objects can also be loaded by a string identifier/name
        json : string
            if already computed a JSON Serialized string of the object
        all_names : dict
            same as cache but for names
        simplifier : util.StorableObjectJSON
            an instance of a JSON Serializer
        identifier : str
            name of the netCDF variable that contains the string to be identified by


        """
        self.storage = storage
        self.content_class = obj
        self.idx_dimension = obj.__name__.lower()
        self.db = obj.__name__.lower()
        self.cache = dict()
        self.named = named
        self.json = json
        self.all_names = None
        self.simplifier = StorableObjectJSON(storage)
        self._names_loaded = False
        if identifier is not None:
            self.identifier = self.idx_dimension + '_' + identifier
        else:
            self.identifier = None

        if dimension_units is not None:
            self.dimension_units = dimension_units
        else:
            self.dimension_units = {}

    @property
    def units(self):
        """
        Return the units dictionary used in the attached storage

        Returns
        -------
        units : dict of {str : simtk.unit.Unit }
            representing a dict of string representing a dimension ('length', 'velocity', 'energy')
            pointing to the simtk.unit.Unit to be used

        """
        return self.storage.units

    def register(self):
        self.storage._storages[self.content_class] = self
        self.storage._storages[self.content_class.__name__] = self
        self.storage._storages[self.content_class.__name__.lower()] = self

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

    def find_by_identifier(self, needle):
        if self.all_names is None:
            self.all_names = { s : idx for idx,s in enumerate(self.storage.variables[self.identifier][:]) }

        if needle in self.all_names:
                return self.all_names[needle]

        return None

    def update_name_cache(self):
        if self.named:
            for idx, name in enumerate(self.variables[self.db + "_name"][:]):
                self.cache[name] = idx

            self._names_loaded = True

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

    def iterator(this, iter_range = None):
        class ObjectIterator:
            def __init__(self):
                self.storage = this
                if iter_range is None:
                    self.idx = 0
                    self.end = self.storage.count()
                else:
                    self.idx = iter_range.start
                    self.end = iter_range.stop

            def __iter__(self):
                return self

            def next(self):
                if self.idx < self.storage.count():
                    obj = self.storage.load(self.idx)
                    self.idx += 1
                    return obj
                else:
                    raise StopIteration()

        return ObjectIterator()

    @loadcache
    def load(self, idx, lazy=True):
        '''
        Returns an object from the storage. Needs to be implented from the specific storage class.
        '''

        # Create object first to break any unwanted recursion in loading
        obj = StoredObject()
        setattr(obj, 'idx', {self.storage : idx})
        if lazy:
            # if lazy construct a function that will update the content. This will be loaded, once the object is accessed
            def loader():
                return self.load_json(self.idx_dimension + '_json', idx, obj)

            setattr(obj, '_loader', loader)
            return obj
        else:
            return self.load_json(self.idx_dimension + '_json', idx, obj)

    @saveidentifiable
    @savecache
    def save(self, obj, idx=None):
        '''
        Returns an object from the storage. Needs to be implented from the specific storage class.
        '''

        if self.named and hasattr(obj, 'name'):
            self.storage.variables[self.db + '_name'][idx] = obj.name

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
        return int(len(self.storage.dimensions[self.idx_dimension]))

    def free(self):
        '''
        Return the number of the next free index

        Returns
        -------
        index : int
            the number of the next free index in the storage. Used to store a new object.
        '''
        return self.count()

    def _init(self, units=None):
        """
        Initialize the associated storage to allow for object storage. Mainly creates an index dimension with the name of the object.

        """
        # define dimensions used for the specific object
        self.storage.createDimension(self.idx_dimension, 0)
        if self.named:
            self.init_variable(self.db + "_name", 'str', description='A short descriptive name for convenience')
        if self.json:
            self.init_variable(self.db + "_json", 'str', description='A json serialized version of the object')

    def var(self, name):
        return '_'.join([self.db, name])

    def begin(self, dimension, idx):
        return int(self.storage.variables[dimension + '_dim_begin'][int(idx)])

    def length(self, dimension, idx):
        return int(self.storage.variables[dimension + '_dim_length'][int(idx)])

    def set_slice(self, dimension, idx, begin, length):
        self.storage.variables[dimension + '_dim_begin'][idx] = begin
        self.storage.variables[dimension + '_dim_length'][idx] = length

    def get_slice(self, dimension, idx):
        begin = int(self.storage.variables[dimension + '_dim_begin'][int(idx)])
        length = int(self.storage.variables[dimension + '_dim_length'][int(idx)])
        return slice(begin, begin+length)

    slice = get_slice

    def free_begin(self, name):
        length = int(len(self.storage.dimensions[name]))
        return length

#=============================================================================================
# INITIALISATION UTILITY FUNCTIONS
#=============================================================================================

    def init_dimension(self, name, length = 0):
        if name not in self.storage.dimensions:
            self.storage.createDimension(name, length)

        self.storage.sync()

    def init_mixed(self, name):
        # entry position is given by a start and a length
        self.init_variable(name + '_dim_begin', 'index', name)
        self.init_variable(name + '_dim_length', 'length', name)

        self.storage.sync()

    def init_variable(self, name, var_type, dimensions = None, units=None, description=None, variable_length=False):
        '''
        Create a new variable in the netcdf storage. This is just a helper function to structure the code better.

        Paramters
        =========
        name : str
            The name of the variable to be created
        var_type : str
            The string representing the type of the data stored in the variable.
            Either the netcdf types can be used directly and
            strings representing the python native types are translated to appropriate netcdf types.
        dimensions : str or tuple of str
            A tuple representing the dimensions used for the netcdf variable.
            If not specified then the default dimension of the storage is used.
        units : str
            A string representing the units used if the var_type is `float` the units is set to `none`
        description : str
            A string describing the variable in a readable form.
        variable_length : bool
            If true the variable is treated as a variable length (list) of the given type. A built-in example
            for this type is a string which is a variable length of char. This make using all the
            mixed stuff superfluous
        '''

        ncfile = self.storage


        if dimensions is None:
            dimensions = self.db

        nc_type = var_type
        if var_type == 'float':
            nc_type = 'f4'   # 32-bit float
        elif var_type == 'int':
            nc_type = np.uint32   # 32-bit signed integer
        elif var_type == 'index':
            nc_type = np.uint32   # 32-bit signed integer / for indices / -1 indicates no index (None)
        elif var_type == 'length':
            nc_type = np.uint32   # 32-bit signed integer / for indices / -1 indicated no length specified (None)
        elif var_type == 'bool':
            nc_type = np.int8   # 8-bit signed integer for boolean
        elif var_type == 'str':
            nc_type = 'str'

        if variable_length:
            vlen_t = self.storage.createVLType(nc_type, name + '_vlen')
            ncvar = ncfile.createVariable(name, vlen_t, dimensions)
        else:
            ncvar = ncfile.createVariable(name, nc_type, dimensions)

        if var_type == 'float' or units is not None:

            unit_instance = u.Unit({})
            symbol = 'none'

            if isinstance(units, u.Unit):
                unit_instance = units
                symbol = unit_instance.get_symbol()
            elif isinstance(units, u.BaseUnit):
                unit_instance = u.Unit({units : 1.0})
                symbol = unit_instance.get_symbol()
            elif type(units) is str and hasattr(u, units):
                unit_instance = getattr(u, units)
                symbol = unit_instance.get_symbol()
            elif type(units) is str and units is not None:
                symbol = units

            json_unit = self.simplifier.unit_to_json(unit_instance)

            # store the unit in the dict inside the Storage object for fast access
            self.storage.units[name] = unit_instance

            # Define units for a float variable
            setattr(ncvar,      'unit_simtk', json_unit)
            setattr(ncvar,      'unit', symbol)

        if description is not None:
            # Define long (human-readable) names for variables.
            setattr(ncvar,    "long_str", description)

        self.storage.sync()


    def init_objectdict(self, name, obj_cls, value_type):
        self.init_dimension(name)

        self.init_variable(name + '_value', value_type, (name, obj_cls))
        self.init_variable(name + '_idx', 'index', (name, obj_cls))
        self.init_variable(name + '_length', 'length', (name))

#=============================================================================================
# LOAD / SAVE UTILITY FUNCTIONS
#=============================================================================================

    def load_variable(self, name, idx):
        return self.storage.variables[name][idx]

    def save_variable(self, name, idx, value):
        self.storage.variables[name][idx] = value

    def load_objectdict(self, name, idx, key_type, value_type):
        length = self.storage.variables[name + '_length'][idx]
        values = self.get_list_as_type(name + '_value', idx, 0, length, value_type)
        keys = self.get_list_as_type(name + '_idx', idx, 0, length, key_type)

        data = dict(zip(keys, values))
        return data

    def save_objectdict(self, name, idx, data, key_type, value_type):
        self.set_list_as_type(name + '_idx', idx, 0, data.keys(), key_type)
        self.set_list_as_type(name + '_value', idx, 0, data.items(), value_type)
        self.save_variable(name + '_length', idx, len(data))

    def load_json(self, name, idx, obj = None):
        idx = int(idx)
        if obj is None:
            obj = StoredObject()

        json_string = self.storage.variables[name][idx]
        setattr(obj, 'json', json_string)

        simplified = yaml.load(json_string)
        data = self.simplifier.build(simplified[2])
        for key, value in data.iteritems():
            setattr(obj, key, value)

        setattr(obj, 'cls', simplified[0])

        return obj

    def save_json(self, name, idx, obj):
        if not hasattr(obj,'json'):
            setattr(obj, 'json', self.object_to_json(obj))

        self.storage.variables[name][idx] = obj.json

#=============================================================================================
# CONVERSION UTILITIES
#=============================================================================================

    def object_to_json(self, obj):
        data = obj.__dict__
        cls = obj.__class__.__name__
#        store = obj.cls
        store = ""

        json_string = self.simplifier.to_json([cls, store, data])

        return json_string

    def list_to_numpy(self, data, value_type, allow_empty = True):
        if value_type == 'int':
            values = np.array(data).astype(np.float32)
        elif value_type == 'float':
            values = np.array(data).astype(np.float32)
        elif value_type == 'bool':
            values = np.array(data).astype(np.int8)
        elif value_type == 'index':
            values = np.array(data).astype(np.uint32)
        elif value_type == 'length':
            values = np.array(data).astype(np.uint32)
        else:
            # an object
            values = [-1 if value is None and allow_empty is True else value.idx[self.storage] for value in data]
            values = np.array(values).astype(np.uint32)

        return values.copy()

    def list_from_numpy(self, values, value_type, allow_empty = True):
        if value_type == 'int':
            data = values.tolist()
        elif value_type == 'float':
            data = values.tolist()
        elif value_type == 'bool':
            data = values.tolist()
        elif value_type == 'index':
            data = values.tolist()
        elif value_type == 'length':
            data = values.tolist()
        else:
            # an object
            key_store = getattr(self.storage, value_type)
            data = [key_store.load(obj_idx) if allow_empty is False or obj_idx >= 0 else None for obj_idx in values.tolist()]

        return data

#=============================================================================================
# SETTER / GETTER UTILITY FUNCTIONS
#=============================================================================================

    def get_object(self, name, idx, cls):
        index = self.load_variable(name + '_idx', idx)
        if index < 0:
            return None

        store = getattr(self.storage, cls)
        obj = store.load(index)
        return obj

    def set_object(self, name, idx, obj):
        if obj is not None:
            self.storage.variables[name + '_idx'][idx] = obj.idx[self.storage]
        else:
            self.storage.variables[name + '_idx'][idx] = -1

    def get_list_as_type(self, name, idx, begin, length, value_type):
        storage = self.storage
        values = storage.variables[name][idx, begin:begin+length]

        data = self.list_from_numpy(values, value_type)
        return data

    def set_list_as_type(self, name, idx, begin, data, value_type):
        values = self.list_to_numpy(data, value_type)
        self.storage.variables[name][idx, begin:begin+len(data)] = values

#=============================================================================================
# ORDERPARAMETER UTILITY FUNCTIONS
#=============================================================================================

    @property
    def op_idx(self):
        """
        Returns aa function that returns for an object of this storage the idx
        """
        def idx(obj):
            return obj.idx[self.storage]

        return idx