import yaml
import types

import numpy as np
import openpathsampling as paths
import simtk.unit as u

import logging
logger = logging.getLogger(__name__)
init_log = logging.getLogger('openpathsampling.initialization')

class Query(object):
    """
    Return
    """
    def __init__(self, query_fnc = None, caching = True):
        if caching:
            self.cache = {}
        else:
            self.cache = None
        self.query_fnc = query_fnc

    def __call__(self, store):
        def _query_iterator(self, store):
            if self.cache is not None:
                if store not in self.cache:
                    self.cache[store] = {}
                store_cache = self.cache[store]

                n_object = len(store)

                for idx in range(n_object):
                    if self.cache is not None:
                        if idx in store_cache:
                            if store_cache[idx]:
                                yield store[idx]
                        else:
                            obj = store.load(idx)
                            result = self.query_fnc(obj)
                            self.cache[store][idx] = result
                            if result:
                                yield obj
                    else:
                        obj = store.load(idx)
                        result = self.query_fnc(obj)
                        if result:
                            yield obj

        return _query_iterator(self, store)

class ObjectStore(object):
    """
    Base Class for storing complex objects in a netCDF4 file. It holds a
    reference to the store file.
    """

    def __init__(self, storage, content_class, has_uid=False, json=True,
                 dimension_units=None, enable_caching=True, load_partial=False,
                 nestable=False, has_name=False):
        """

        Parameters
        ----------
        storage
        content_class
        has_uid
        json
        dimension_units
        enable_caching : bool
            if this is set to `True` caching is used to quickly access
            previously loaded objects (default)
        load_partial : bool
            if this is set to `True` the storage allows support for partial
            delayed loading of member variables. This is useful for larger
            objects that might only be required in particular circumstances.
            (default is `False`)
        nestable : bool
            if true this marks the content_class to be saved as nested dict
            objects and not a pointing to saved objects. So the saved complex
            object is only stored once and not split into several objects that
            are referenced by each other in a tree-like fashion


        Attributes
        ----------

        storage : Storage
            the reference the Storage object where all data is stored
        content_class : class
            a reference to the class type to be stored using this Storage
        has_uid : bool
            if `True` objects can also be loaded by a string identifier/name
        json : string
            if already computed a JSON Serialized string of the object
        dimension_units : dict of {str : simtk.unit.Unit } or None
            representing a dict of string representing a dimension
            ('length', 'velocity', 'energy') pointing to
            the simtk.unit.Unit to be used. If not None overrides the standard
            units used in the storage
        simplifier : util.StorableObjectJSON
            an instance of a JSON Serializer
        identifier : str
            name of the netCDF variable that contains the string to be
            identified by. So far this is `name`
        cache : dict (int or str : object)
            a dictionary that holds references to all stored elements by index
            or string for named objects. This is only used for cached access
            is enable_caching is True (default)

        Notes
        -----
        The class that takes care of storing data in a file is called a Storage,
        so the netCDF subclassed Storage is a storage. The classes that know how
        to load and save an object from the storage are called stores,
        like ObjectStore, SampleStore, etc...

        """
        self.storage = storage
        self.content_class = content_class
        self.idx_dimension = content_class.__name__.lower()
        self.db = content_class.__name__.lower()
        self.cache = dict()
        self.has_uid = has_uid
        self.has_name = has_name
        self.json = json
        self.simplifier = paths.storage.StorableObjectJSON(storage)
        self.identifier = self.db + '_uid'
        self._free = set()
        self._cached_all = False
        self._names_loaded = False

        if dimension_units is not None:
            self.dimension_units = dimension_units
        else:
            self.dimension_units = self.storage.dimension_units

        # First, apply standard decorator for loading and saving
        # this handles all the setting and getting of .idx and is
        # always necessary!

        if load_partial:
            # this allows the class to load members only if needed
            # adds a different __getattr__ to the content class
            # makes only sense if not already lazy loading
            # it uses load_constructor instead to create an empty object
            # and then each class can attach delayed loaders to load
            # when necessary, fall back is of course the normal load function

            if hasattr(self, 'load_empty'):
                cls = self.content_class

                def _getattr(this, item):
                    if item == '_idx':
                        return this.__dict__['idx']

                    if hasattr(cls, '_delayed_loading'):
                        if item in dir(cls):
                            return object.__getattribute__(this, item)

                        if item in cls._delayed_loading:
                            _loader = cls._delayed_loading[item]
#                            print 'from', repr(self.storage), id(self), 'and not', repr(this), 'load', item
                            _loader(this)
                        else:
                            raise KeyError(item)

                    return this.__dict__[item]

                setattr(cls, '__getattr__', _getattr)

                _load = self.load
                self.load = types.MethodType(loadpartial(_load), self)

        _save = self.save
        self.save = types.MethodType(saveidx(_save), self)

        _load = self.load
        self.load = types.MethodType(loadidx(_load), self)

        if enable_caching:
            # wrap load/save to make this work. I use MethodType here to bind the
            # wrapped function to this instance. An alternative would be to
            # add the wrapper to the class itself, which would mean that all
            # instances have the same setting and a change would be present in all
            # instances. E.g., first instance has caching and the second not
            # when the second instance is created the change in the class would
            # also disable caching in the first instance. The present way with
            # bound methods is more flexible
            # Should be not really important, since there will be mostly only one
            # storage, but this way it is cleaner

            _save = self.save
            self.save = types.MethodType(savecache(_save), self)

            _load = self.load
            self.load = types.MethodType(loadcache(_load), self)

        self._register_with_storage(nestable=nestable)

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return "%s(content_class=%s, variable_prefix=%s)" % (
            self.__class__.__name__, self.content_class, self.db)


    def idx(self, obj):
        """
        Return the index in this store for a given object

        Parameters
        ----------
        obj : object
            the object that can be stored in this store for which its index is
            to be returned

        Returns
        -------
        int or None
            The integer index of the given object or None if it is not stored yet
        """
        if hasattr(obj, 'idx'):
            if self.storage in obj.idx:
                return obj.idx[self.storage]

        return None

    @property
    def units(self):
        """
        Return the units dictionary used in the attached storage

        Returns
        -------
        units : dict of {str : simtk.unit.Unit }
            representing a dict of string representing a dimension
            ('length', 'velocity', 'energy')
            pointing to the simtk.unit.Unit to be used

        """
        return self.storage.units

    def _register_with_storage(self, nestable = False):
        """
        Register the store with the associated storage. This means the

        Parameters
        ----------
        nestable : bool
            if true this marks the content_class to be saved as nested dict
            objects and not a pointing to saved objects. So the saved complex
            object is only stored once and not split into several objects that
            are referenced by each other in a tree-like fashion

        Notes
        -----
        """

        self.storage._storages[self.content_class] = self
        self.storage._storages[self.content_class.__name__] = self
        self.storage._storages[self.content_class.__name__.lower()] = self

        # Add here the logic to change the actual class and add the decorator
        # this is the same as add the  decorator (without default_storage,
        # I removed this since it is not used anyway)
        # all it does is make sure that there is a .idx property and the base_
        # cls is known
        self.content_class.base_cls_name = self.content_class.__name__
        self.content_class.base_cls = self.content_class

        # add a property idx that keeps the storage reference
        def _idx(this):
            if not hasattr(this, '_idx'):
                this._idx = dict()

            return this._idx

        self.content_class.idx = property(_idx)

        if self.has_uid:
            def _uid_get(this):
                if not hasattr(this, '_uid'):
                    this._uid = ''

                return this._uid

            def _uid_set(this, uid):
                this._uid = uid

            self.content_class.uid = property(_uid_get, _uid_set)

        def _save(this, storage):
            if storage is not None:
                storage.save(this)

        if nestable:
            self.content_class.nestable = True

        self.content_class.save = _save

        if not hasattr(self.content_class, 'cls'):
            def _cls(this):
                return this.__class__.__name__

            self.content_class.cls = property(_cls)

        # register as a base_class for storable objects
        self.storage.links.append(self)


    def set_variable_partial_loading(self, variable, loader):
        cls = self.content_class
        if not hasattr(cls, '_delayed_loading'):
            cls._delayed_loading = dict()

        cls._delayed_loading[variable] = loader

    def idx_by_name(self, needle):
        """
        Return the index for the (first) object with a given name from the store

        Parameters
        ----------
        needle : str
            The name of the object to be found in the storage

        Returns
        -------
        int or None
            The index of the first found object. If the name is not present,
            None is returned

        Notes
        -----
        Can only be applied to named storages.
        """
        if self.has_uid:
            # if we need a cache we might find the index in there
            if needle in self.cache:
                if type(self.cache[needle]) is int:
                    return self.cache[needle]
                else:
                    return self.cache[needle].idx[self.storage]

            # otherwise search the storage for the name
            found_idx = [ idx for idx,s in enumerate(self.storage.variables[
                self.identifier][:]) if s == needle
            ]

            if len(found_idx) > 0:
                    return found_idx[0]

            return None
        else:
            raise ValueError('Cannot search for name (str) in non-named objects')

    def update_name_cache(self):
        """
        Update the internal cache with all stored names in the store.
        This allows to load by name for named objects
        """
        if self.has_name:
            if not self._names_loaded:
                for idx, name in enumerate(self.storage.variables[self.db + "_name"][:]):
                    self._update_name_in_cache(name, idx)

                self._names_loaded = True

    def _update_name_in_cache(self, name, idx):
        if name != '':
            if name not in self.cache:
                self.cache[name] = [idx]
            else:
                if idx not in self.cache[name]:
                    self.cache[name].append(idx)

    def find(self, name):
        """
        Return all objects with a given name

        Parameters
        ----------
        name : str
            the name to be searched for

        Returns
        -------
        list of objects
            a list of found objects, can be empty [] if no objects with
            that name exist

        """
        if self.has_name:
            if name in self.cache:
                self.update_name_cache()

            return self[self.cache[name]]

        return []

    def find_indices(self, name):
        """
        Return indices for all objects with a given name

        Parameters
        ----------
        name : str
            the name to be searched for

        Returns
        -------
        list of int
            a list of indices in the storage for all found objects,
            can be empty [] if no objects with that name exist

        """
        if self.has_name:
            if name in self.cache:
                self.update_name_cache()

            return self.cache[name]

        return []


    def find_first(self, name):
        """
        Return first object with a given name

        Parameters
        ----------
        name : str
            the name to be searched for

        Returns
        -------
        object of None
            the first found object, can be None if no object with the given
            name exists

        """
        if self.has_name:
            if name in self.cache:
                self.update_name_cache()

            if len(self.cache[name]) > 0:
                return self[self.cache[name][0]]

        return None


    def __iter__(self):
        """
        Add iteration over all elements in the storage
        """
        return self.iterator()

    def __len__(self):
        """
        Return the number of stored objects

        Returns
        -------
        int
            number of stored objects

        Notes
        -----
        Equal to `store.count()`
        """
        return self.count()

    def iterator(this, iter_range = None):
        """
        Return an iterator over all objects in the storage

        Parameters
        ----------
        iter_range : slice or None
            if this is not `None` it confines the iterator to objects specified
            in the slice

        Returns
        -------
        Iterator()
            The iterator that iterates the objects in the store

        """
        class ObjectIterator:
            def __init__(self):
                self.storage = this
                self.iter_range = iter_range
                if iter_range is None:
                    self.idx = 0
                    self.end = self.storage.count()
                else:
                    self.idx = iter_range.start
                    self.end = iter_range.stop

            def __iter__(self):
                return self

            def next(self):
                if self.idx < self.end:
                    obj = self.storage.load(self.idx)
                    if self.iter_range is not None and self.iter_range.step is not None:
                        self.idx += self.iter_range.step
                    else:
                        self.idx += 1
                    return obj
                else:
                    raise StopIteration()

        return ObjectIterator()

    def __getitem__(self, item):
        """
        Enable numpy style selection of object in the store
        """
        try:
            if type(item) is int or type(item) is str:
                return self.load(item)
            elif type(item) is slice:
                return [self.load(idx) for idx in range(*item.indices(len(self)))]
            elif type(item) is list:
                return [self.load(idx) for idx in item]
            elif item is Ellipsis:
                return self.iterator()
        except KeyError:
            return None

    def load(self, idx):
        '''
        Returns an object from the storage. Needs to be implemented from
        the specific storage class.

        Parameters
        ----------
        idx : int or str
            either the integer index of the object to be loaded or a string
            (name) for named objects. This will always return the first object
            found with the specified name.

        Returns
        -------
        object
            the loaded object
        '''

        return self.load_json(self.idx_dimension + '_json', idx)

    def clear_cache(self):
        """Clear the cache and force reloading

        """

        self.cache = dict()
        self._cached_all = False

    def cache_all(self):
        """Load all samples as fast as possible into the cache

        """
        if not self._cached_all:
            idxs = range(len(self))
            jsons = self.storage.variables[self.idx_dimension + '_json'][:]

            [ self.add_single_to_cache(i,j) for i,j in zip(
                idxs,
                jsons) ]

            self._cached_all = True

    def add_single_to_cache(self, idx, json):
        """
        Add a single object to cache by json
        """

        simplified = yaml.load(json)
        obj = self.simplifier.build(simplified)

        obj.json = json
        obj.idx[self.storage] = idx

        self.cache[idx] = obj

        return obj


    def save(self, obj, idx=None):
        """
        Saves an object the storage.

        Parameters
        ----------
        obj : object
            the object to be stored
        idx : int or string or `None`
            the index to be used for storing. This is highly discouraged since
            it changes an immutable object (at least in the storage). It is
            better to store also the new object and just ignore the
            previously stored one.

        """

        if self.has_uid and hasattr(obj, '_uid'):
            self.storage.variables[self.identifier][idx] = obj._uid

        self.save_json(self.idx_dimension + '_json', idx, obj)

    def get_uid(self, idx):
        """
        Return the name of and object with given integer index

        Parameters
        ----------
        idx : int
            the integer index of the object whose name is to be returned

        Returns
        -------
        str or None
            Returns the name of the object for named objects. None otherwise.

        """
        if self.has_uid:
            return self.storage.variables[self.identifier][idx]
        else:
            return None

    def get(self, indices):
        """
        Returns a list of objects from the given list of indices

        Parameters
        ----------
        indices : list of int
            the list of integers specifying the object to be returned

        Returns
        -------
        list of objects
            a list of objects stored under the given indices

        """
        return [self.load(idx) for idx in range(0, self.count())[indices]]

    def last(self):
        '''
        Returns the last generated trajectory. Useful to continue a run.

        Returns
        -------
        Trajectoy
            the actual trajectory object
        '''
        return self.load(self.count() - 1)

    def first(self):
        '''
        Returns the last stored object. Useful to continue a run.

        Returns
        -------
        Object
            the actual last stored object
        '''
        return self.load(0)

    def count(self):
        '''
        Return the number of objects in the storage

        Returns
        -------
        number : int
            number of objects in the storage.

        Notes
        -----
        Use len(store) instead
        '''
        return int(len(self.storage.dimensions[self.idx_dimension]))

    def free(self):
        '''
        Return the number of the next free index

        Returns
        -------
        index : int
            the number of the next free index in the storage.
            Used to store a new object.
        '''
        count = self.count()
        self._free = set([ idx for idx in self._free if idx >= count])
        idx = count
        while idx in self._free:
            idx += 1

        return idx

    def reserve_idx(self, idx):
        '''
        Locks an idx as used
        '''
        self._free.add(idx)

    def _init(self, units=None):
        """
        Initialize the associated storage to allow for object storage. Mainly
        creates an index dimension with the name of the object.

        Parameters
        ----------
        units : dict of {str : simtk.unit.Unit} or None
            representing a dict of string representing a dimension
            ('length', 'velocity', 'energy') pointing to
            the simtk.unit.Unit to be used. If not None overrides the standard
            units used in the storage
        """
        # define dimensions used for the specific object
        self.storage.createDimension(self.idx_dimension, 0)
#        if self.has_name:
#            self.init_variable(self.db + "_name", 'str',
#                description='A short descriptive name for convenience',
#                chunksizes=tuple([10240]))

        if self.has_uid:
            self.init_variable(self.db + "_uid", 'str',
                description='A unique identifier',
                chunksizes=tuple([10240]))

        if self.has_name:
            self.init_variable(self.db + "_name", 'str',
                description='A name',
                chunksizes=tuple([10240]))

        if self.json:
            self.init_variable(self.db + "_json", 'str',
                description='A json serialized version of the object',
                chunksizes=tuple([10240]))

#==============================================================================
# INITIALISATION UTILITY FUNCTIONS
#==============================================================================

    def init_dimension(self, name, length = 0):
        """
        Initialize a new dimension in the storage.
        Wraps the netCDF createDimension

        Parameters
        ----------
        name : str
            the name for the new dimension
        length : int
            the number of elements in this dimension. Zero (0) (default) means
            an infinite dimension that extends when more objects are stored

        """
        if name not in self.storage.dimensions:
            self.storage.createDimension(name, length)

#        self.storage.sync()

    @staticmethod
    def _parse_var_type_as_np_type(var_type):
        nc_type = var_type
        if var_type == 'float':
            nc_type = np.float32   # 32-bit float
        elif var_type == 'int':
            nc_type = np.int32   # 32-bit signed integer
        elif var_type == 'index':
            nc_type = np.int32
            # 32-bit signed integer / for indices / -1 : no index (None)
        elif var_type == 'length':
            nc_type = np.int32
            # 32-bit signed integer / for indices / -1 : no length specified (None)
        elif var_type == 'bool':
            nc_type = np.uint8   # 8-bit signed integer for boolean
        elif var_type == 'str':
            nc_type = 'str'

        types = {
            'float' : np.float32,
            'int' : np.int32,
            'index' : np.int32,
            'length' : np.int32,
            'bool' : np.uint8,
            'str' : 'str',
            float : np.float32,
            int : np.int32,
            bool : np.int8,
            str : 'str',
            np.uint8 : np.uint8,
            np.uint16 : np.uint16,
            np.uint32 : np.uint32,
            np.uint64 : np.uint64,
            np.int8 : np.int8,
            np.int16 : np.int16,
            np.int32 : np.int32,
            np.int64 : np.int64,
            np.float32 : np.float32,
            np.float64 : np.float64
        }

        return types[var_type]

    def init_variable(self, name, var_type, dimensions = None, units=None,
                      description=None, variable_length=False, chunksizes=None):
        '''
        Create a new variable in the netCDF storage. This is just a helper
        function to structure the code better.

        Parameters
        ==========
        name : str
            The name of the variable to be created
        var_type : str
            The string representing the type of the data stored in the variable.
            Either the netCDF types can be used directly and
            strings representing the python native types are translated to
            appropriate netCDF types.
        dimensions : str or tuple of str
            A tuple representing the dimensions used for the netcdf variable.
            If not specified then the default dimension of the storage is used.
        units : str
            A string representing the units used if the var_type is `float`
            the units is set to `none`
        description : str
            A string describing the variable in a readable form.
        variable_length : bool
            If true the variable is treated as a variable length (list) of the
            given type. A built-in example for this type is a string which is
            a variable length of char. This make using all the mixed
            stuff superfluous
        chunksizes : tuple of int
            A tuple of ints per number of dimensions. This specifies in what
            block sizes a variable is stored. Usually for object related stuff
            we want to store everything of one object at once so this is often
            (1, ..., ...)
        '''

        ncfile = self.storage

        if dimensions is None:
            dimensions = self.db

        nc_type = self._parse_var_type_as_np_type(var_type)

        if variable_length:
            vlen_t = ncfile.createVLType(nc_type, name + '_vlen')
            ncvar = ncfile.createVariable(name, vlen_t, dimensions,
                                          zlib=False, chunksizes=chunksizes)
        else:
            ncvar = ncfile.createVariable(name, nc_type, dimensions,
                                          zlib=False, chunksizes=chunksizes)

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

            # store the unit in the dict inside the Storage object
            self.storage.units[name] = unit_instance

            # Define units for a float variable
            setattr(ncvar,      'unit_simtk', json_unit)
            setattr(ncvar,      'unit', symbol)

        if description is not None:
            # Define long (human-readable) names for variables.
            setattr(ncvar,    "long_str", description)

#        self.storage.sync()

#==============================================================================
# LOAD / SAVE UTILITY FUNCTIONS
#==============================================================================

    def load_variable(self, name, idx):
        """
        Wrapper for netCDF storage.variables[name][idx] property

        Parameters
        ----------
        name : str
            The name of the variable
        idx : int, slice, list of int, etc...
            An index specification as in netCDF4

        Returns
        -------
        numpy.ndarray
            The data stored in the netCDF variable

        """
        return self.storage.variables[name][idx]

    def save_variable(self, name, idx, value):
        """
        Wrapper for netCDF storage.variables[name][idx] property

        Parameters
        ----------
        name : str
            The name of the variable
        idx : int, slice, list of int, etc...
            An index specification as in netCDF4
        value : numpy.ndarray
            The array to be stored in the variable

        """
        self.storage.variables[name][idx] = value

    def load_json(self, name, idx):
        """
        Load an object from the associated storage using json

        Parameters
        ----------
        name : str
            the name of the variable in the netCDF storage
        idx : int
            the integer index in the variable

        Returns
        -------
        object
            the loaded object

        """
        # TODO: Add logging here
        idx = int(idx)

        json_string = self.storage.variables[name][idx]

        simplified = yaml.load(json_string)
        obj = self.simplifier.build(simplified)
        setattr(obj, 'json', json_string)

        return obj

    def save_json(self, name, idx, obj):
        """
        Save an object as a json string in a variable in the referenced storage

        Parameters
        ----------
        name : str
            the name of the variable in the netCDF storage
        idx : int
            the integer index in the variable
        obj : object
            the object to be stored as JSON

        """
        if not hasattr(obj,'json'):
            setattr(obj, 'json', self.object_to_json(obj))

        self.storage.variables[name][idx] = obj.json


#==============================================================================
# CONVERSION UTILITIES
#==============================================================================

    def object_to_json(self, obj):
        """
        Convert a given object to a json string using the simplifier

        Parameters
        ----------
        obj : the object to be converted

        Returns
        -------
        str
            the JSON string
        """
        json_string = self.simplifier.to_json_object(obj, obj.base_cls_name)

        return json_string

    def list_to_numpy(self, data, value_type, allow_empty = True):
        """
        Return a numpy list from a python list in a given format

        Parameters
        ----------
        data : list
            the list to be converted
        value_type : str
            the type of the input list elements. If this is an object type it
            will be saved and the returned index is stored in an numpy
            integer array
        allow_empty : bool
            if set to `True` None will be stored as the integer -1

        Returns
        -------
        numpy.ndarray
            the converted numpy array
        """
        if value_type == 'int':
            values = np.array(data).astype(np.float32)
        elif value_type == 'float':
            values = np.array(data).astype(np.float32)
        elif value_type == 'bool':
            values = np.array(data).astype(np.uint8)
        elif value_type == 'index':
            values = np.array(data).astype(np.int32)
        elif value_type == 'length':
            values = np.array(data).astype(np.int32)
        else:
            # an object
            values = [-1 if value is None and allow_empty is True
                      else value.idx[self.storage] for value in data]
            values = np.array(values).astype(np.int32)

        return values.copy()

    def list_from_numpy(self, values, value_type, allow_empty = True):
        """
        Return a python list from a numpy array in a given format

        Parameters
        ----------
        values : numpy.ndarray
            the numpy array to be converted
        value_type : str
            the type of the output list elements. If this is a object type it
            will be loaded using the numpy array content as the index
        allow_empty : bool
            if set to `True` then loaded objects will only be loaded if the
            index is not negative. Otherwise the load function will always
            be called

        Returns
        -------
        list
            the converted list
        """
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
            data = [key_store.load(obj_idx) if allow_empty is False
                    or obj_idx >= 0 else None for obj_idx in values.tolist()]

        return data

#==============================================================================
# SETTER / GETTER UTILITY FUNCTIONS
#==============================================================================

    # TODO: This might go tho storage.py
    def load_object(self, name, idx, store):
        """
        Load an object from the storage

        Parameters
        ----------
        name : str
            name of the variable to be used
        index : int
            index in the storage
        cls : cls
            type of the object to be loaded. Determines the store to be used

        Returns
        -------
        object
            the loaded object
        """
        index = self.load_variable(name + '_idx', idx)
        if index < 0:
            return None

        obj = store.load(index)
        return obj

    def save_object(self, name, idx, obj):
        """
        Store an object in the storage

        Parameters
        ----------
        name : str
            name of the variable to be used
        index : int
            index in the storage
        obj : object
            the object to be stored

        """
        storage = self.storage

        if obj is not None:
            storage.save(obj)
            storage.variables[name + '_idx'][idx] = obj.idx[storage]
        else:
            storage.variables[name + '_idx'][idx] = -1

#==============================================================================
# COLLECTIVE VARIABLE UTILITY FUNCTIONS
#==============================================================================

    @property
    def op_idx(self):
        """
        Returns a function that returns for an object of this storage the idx.
        This can be used to construct order parameters the return the index
        in this storage. Useful for visualization

        Returns
        -------
        function
            the function that reports the index in this store
        """
        def idx(obj):
            return obj.idx[self.storage]

        return idx

#=============================================================================
# LOAD/SAVE DECORATORS FOR PARTIAL LOADING OF ATTRIBUTES
#=============================================================================

def loadpartial(func, constructor=None):
    """
    Decorator for load functions that add the basic handling for partial loading
    """

    def inner(self, idx, *args, **kwargs):
        if constructor is None:
            new_func = getattr(self, 'load_empty')
        else:
            new_func = getattr(self, constructor)

        return_obj = new_func(idx, *args, **kwargs)
        # this tells the obj where it was loaded from
        return_obj._origin = self.storage
        return return_obj

    return inner


#=============================================================================
# LOAD/SAVE DECORATORS FOR CACHE HANDLING
#=============================================================================

def loadcache(func):
    """
    Decorator for load functions that add the basic cache handling
    """
    def inner(self, idx, *args, **kwargs):
        if type(idx) is not str and idx < 0:
            return None

        n_idx = idx

        # if it is in the cache, return it, otherwise not :)
        if idx in self.cache:

            cc = self.cache[idx]
            if type(cc) is int:
                logger.debug('Found IDX #' + str(idx) + ' in cache under position #' + str(cc))

                # here the cached value is actually only the index
                # so it still needs to be loaded with the given index
                # this happens when we want to load by name (str)
                # and we need to actually load it
                n_idx = cc
            elif type(cc) is list:
                logger.debug('Found IDX #' + str(idx) + ' in cache under positions #' + str(cc) + '. Loading first.')
                n_idx = cc[0]

            else:
                logger.debug('Found IDX #' + str(idx) + ' in cache. Not loading!')

                # we have a real object (hopefully) and just return from cache
                n_idx = idx

            if n_idx in self.cache:
                # return from cache
                return self.cache[n_idx]


        elif type(idx) is str:
            # we want to load by name and it was not in cache.
            if self.has_name:
                # since it is not found in the cache before. Refresh the cache
                self.update_name_cache()

                # and give it another shot
                if idx in self.cache:
                    if len(self.cache[idx]) > 1:
                        logger.debug('Found name "%s" multiple (%d) times in storage! Loading first!' % (idx, len(self.cache[idx])))
                        n_idx = self.cache[idx][0]
                    else:
                        n_idx = self.cache[idx][0]
                else:
                    raise ValueError('str "' + idx + '" not found in storage')
            else:
                raise ValueError('str "' + idx + '" as indices are only allowed in named storage')

        # ATTENTION HERE!
        # Note that the wrapped function no self as first parameter. This is because we are wrapping a bound
        # method in an instance and this one is still bound - luckily - to the same 'self'. In a class decorator when wrapping
        # the class method directly it is not bound yet and so we need to include the self! Took me some time to
        # understand and figure that out
        obj = func(n_idx, *args, **kwargs)
        if obj is not None:
            # update cache there might have been a change due to naming
            self.cache[obj.idx[self.storage]] = obj

            # finally store the name of a named object in cache
            if self.has_name and obj._name != '':
                self._update_name_in_cache(obj._name, n_idx)

        return obj
    return inner

# the default decorator for save functions to enable caching
def savecache(func):
    """
    Decorator for save functions that add the basic cache handling
    """
    def inner(self, obj, idx = None, *args, **kwargs):
        # call the normal storage
        func(obj, idx, *args, **kwargs)
        idx = obj.idx[self.storage]

        # store the name in the cache
        self.cache[idx] = obj
        if self.has_name and obj._name != '':
            # and also the name, if it has one so we can load by
            # name afterwards from cache
                self._update_name_in_cache(obj._name, idx)

    return inner

#=============================================================================
# LOAD/SAVE DECORATORS FOR .idx HANDLING
#=============================================================================

def loadidx(func):
    """
    Decorator for load functions that add the basic indexing handling
    """
    def inner(self, idx, *args, **kwargs):
        if type(idx) is not str and int(idx) < 0:
            return None

        n_idx = idx

        if type(idx) is str:
            # we want to load by name and it was not in cache
            if self.has_name:
                raise ValueError('Load by name without caching is not supported')
#                n_idx = self.load_by_name(idx)
            else:
                # load by name only in named storages
                raise ValueError('Load by name (str) is only supported in named storages')
                pass

        # turn into python int if it was a numpy int (in some rare cases!)
        n_idx = int(n_idx)

        # ATTENTION HERE!
        # Note that the wrapped function ho self as first parameter. This is because we are wrapping a bound
        # method in an instance and this one is still bound - luckily - to the same 'self'. In a class decorator when wrapping
        # the class method directly it is not bound yet and so we need to include the self! Took me some time to
        # understand and figure that out
        logger.debug('Calling load object of type ' + self.content_class.__name__ + ' and IDX #' + str(idx))
        if n_idx >= len(self):
            logger.warning('Trying to load from IDX #' + str(n_idx) + ' > number of object ' + str(len(self)))
            return None
        elif n_idx < 0:
            logger.warning('Trying to load negative IDX #' + str(n_idx) + ' < 0')
            return None
        else:
            obj = func(n_idx, *args, **kwargs)

        if not hasattr(obj, 'idx'):
            obj.idx = dict()

        obj.idx[self.storage] = n_idx

        if self.has_uid:
            if not hasattr(obj, '_uid'):
                # get the name of the object
                setattr(obj, '_uid', self.get_uid(idx))

        if self.has_name and hasattr(obj, '_name'):
            setattr(obj, '_name',
                    self.storage.variables[self.db + '_name'][idx])
            # make sure that you cannot change the name of loaded objects
            obj.fix_name()

        return obj

    return inner

def saveidx(func):
    """
    Decorator for save functions that add the basic indexing handling
    """
    def inner(self, obj, idx = None, *args, **kwargs):

        storage = self.storage
        if idx is None:
            if storage in obj.idx:
                # has been saved so quit and do nothing
                return None
            else:
                idx = self.free()
        else:
            if type(idx) is str:
                # Not yet supported
                raise ValueError('Saving by name not yet supported')
            else:
                idx = int(idx)

        obj.idx[storage] = idx

        # make sure in nested saving that an IDX is not used twice!
        self.reserve_idx(idx)
        logger.debug('Saving ' + str(type(obj)) + ' using IDX #' + str(idx))
        func(obj, idx, *args, **kwargs)

        if self.has_uid and hasattr(obj, '_uid') and obj._uid != '':
            self.storage.variables[self.identifier][idx] = obj._uid

        if self.has_name and hasattr(obj, '_name'):
            #logger.debug('Object ' + str(type(obj)) + ' with IDX #' + str(idx))
            #logger.debug(repr(obj))
            #logger.debug("Cleaning up name; currently: " + str(obj._name))
            if obj._name is None:
                # this should not happen!
                logger.debug("Nameable object has not been initialized correctly. Has None in _name")
                raise AttributeError('_name needs to be a string for nameable objects.')

            obj.fix_name()

            self.storage.variables[self.db + '_name'][idx] = obj._name

    return inner

# CREATE EASY UPDATE WRAPPER

def func_update_object(attribute, db, variable, store):
    """
    Create a delayed loading function for stores

    Parameters
    ----------
    attribute : string
        name of the attribute of the object to be updated. E.g. for sample.mover this is 'mover'
    db : string
        the storage prefix where the object are stored in the file. E.g. for samples this is 'sample'
    variable : string
        the name of the variable in the storage. this is often the same as the attribute
    store : string
        the name of the store. E.g. 'trajectories'

    Returns
    -------
    function
        the function that is used for updating
    """
    def updater(obj):
        storage = obj._origin

        idx = obj.idx[storage]
        obj_idx = int(storage.variables[db + '_' + variable + '_idx'][idx])

        setattr(obj, attribute, getattr(storage, store).load(obj_idx))

    return staticmethod(updater)
