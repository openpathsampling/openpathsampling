import types
import logging

import yaml

from cache import MaxCache, Cache, NoCache, WeakLRUCache
from objproxy import LoaderProxy

import weakref

logger = logging.getLogger(__name__)
init_log = logging.getLogger('openpathsampling.initialization')


class ObjectStore(object):
    """
    Base Class for storing complex objects in a netCDF4 file. It holds a
    reference to the store file.
    """

    allowed_types = [
        'int', 'float', 'long', 'str', 'bool'
                                       'numpy.float32', 'numpy.float64',
        'numpy.int8', 'numpy.inf16', 'numpy.int32', 'numpy.int64',
        'numpy.uint8', 'numpy.uinf16', 'numpy.uint32', 'numpy.uint64',
        'index', 'length'
    ]

    class DictDelegator(object):
        def __init__(self, store, dct):
            self.prefix = store.prefix + '_'
            self.dct = dct

        def __getitem__(self, item):
            return self.dct[self.prefix + item]

    def prefix_delegate(self, dct):
        return ObjectStore.DictDelegator(self, dct)

    default_cache = 10000

    def __init__(self, content_class, json=True,
                 caching=None, nestable=False, has_name=False):

        """

        Parameters
        ----------
        storage
        content_class
        json
        dimension_units
        caching : dict-like or bool or int or None
            this is the dict used for caching.
            `True` means to use a python built-in dict which unlimited caching.
            Be careful.
            `False` means no caching at all. If a dict-like object is passed,
            it will be used.
            An integer `n` means to use LRU Caching with maximal n elements and is
            equal to `cache=LRUCache(n)`
            Default (None) is equivalent to `cache=ObjectStore.default_cache`
        nestable : bool
            if true this marks the content_class to be saved as nested dict
            objects and not a pointing to saved objects. So the saved complex
            object is only stored once and not split into several objects that
            are referenced by each other in a tree-like fashion

        Notes
        -----
        Usually you want caching, but limited. Recommended is to use an LRUCache
        with a reasonable number that depends on the typical number of objects to
        cache and their size

        Attributes
        ----------

        storage : Storage
            the reference the Storage object where all data is stored
        content_class : class
            a reference to the class type to be stored using this Storage
        has_name : bool
            if `True` objects can also be loaded by a string name
        json : string
            if already computed a JSON Serialized string of the object
        simplifier : util.StorableObjectJSON
            an instance of a JSON Serializer
        identifier : str
            name of the netCDF variable that contains the string to be
            identified by. So far this is `name`
        cache : dict-like (int or str : object)
            a dictionary that holds references to all stored elements by index
            or string for named objects. This is only used for cached access
            if caching is not `False`

        Notes
        -----
        The class that takes care of storing data in a file is called a Storage,
        so the netCDF subclassed Storage is a storage. The classes that know how
        to load and save an object from the storage are called stores,
        like ObjectStore, SampleStore, etc...

        """

        self._storage = None
        self.content_class = content_class
        self.prefix = None
        self.cache = NoCache()
        self.has_name = has_name
        self.json = json
        self._free = set()
        self._cached_all = False
        self._names_loaded = False
        self.nestable = nestable
        self.name_idx = dict()

        self.variables = dict()
        self.vars = dict()
        self.units = dict()

        self.index = weakref.WeakKeyDictionary()

    def register(self, storage, name):
        self._storage = storage
        self.prefix = name

        self.variables = self.prefix_delegate(self.storage.variables)
        self.units = self.prefix_delegate(self.storage.units)
        self.vars = self.prefix_delegate(self.storage.vars)

    @property
    def storage(self):
        if self._storage is None:
            raise RuntimeError('A store need to be added to a storage to be used!')

        return self._storage

    @property
    def dimension_units(self):
        return self.storage.dimension_units

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return "store.%s[%s]" % (
            self.prefix, self.content_class.__name__)

    @property
    def simplifier(self):
        return self.storage.simplifier

    def set_caching(self, caching):
        if caching is None:
            caching = self.default_cache

        if caching is True:
            caching = MaxCache()
        elif caching is False:
            caching = NoCache()
        elif type(caching) is int:
            caching = WeakLRUCache(caching)

        if isinstance(caching, Cache):
            self.cache = caching.transfer(self.cache)

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
        return self.index.get(obj, None)

    def update_name_cache(self):
        """
        Update the internal cache with all stored names in the store.
        This allows to load by name for named objects
        """
        if self.has_name:
            if not self._names_loaded:
                for idx, name in enumerate(self.storage.variables[self.prefix + "_name"][:]):
                    self._update_name_in_cache(name, idx)

                self._names_loaded = True

    def _update_name_in_cache(self, name, idx):
        if name != '':
            if name not in self.cache:
                self.name_idx[name] = [idx]
            else:
                if idx not in self.cache[name]:
                    self.name_idx[name].append(idx)

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
            if name not in self.name_idx:
                self.update_name_cache()

            return self[self.name_idx[name]]

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
            if name not in self.name_idx:
                self.update_name_cache()

            return self.name_idx[name]

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
            if name not in self.name_idx:
                self.update_name_cache()

            if len(self.name_idx[name]) > 0:
                return self[self.name_idx[name][0]]

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
        return len(self.storage.dimensions[self.prefix])

    def iterator(this, iter_range=None):
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
                    self.end = len(self.storage)
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

    def write(self, variable, idx, obj, attribute=None):
        if attribute is None:
            attribute = variable

        var = self.vars[variable]
        val = getattr(obj, attribute)

        var[int(idx)] = val

        if var.var_type.startswith('lazy'):
            proxy = var.store.proxy(val)
            setattr(obj, attribute, proxy)

    def proxy(self, item):
        if item is None:
            return None

        if type(item) is not int:
            idx = self.index.get(item, None)
        else:
            idx = item

        if idx is None:
            return item
        else:
            return LoaderProxy(self, idx)

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

    def _load(self, idx):
        obj = self.vars['json'][idx]
        return obj

    #        return self.load_json(self.prefix + '_json', idx)

    def clear_cache(self):
        """Clear the cache and force reloading

        """

        self.cache.clear()
        self._cached_all = False

    def cache_all(self):
        """Load all samples as fast as possible into the cache

        """
        if not self._cached_all:
            idxs = range(len(self))
            jsons = self.variables['json'][:]

            [self.add_single_to_cache(i, j) for i, j in zip(
                idxs,
                jsons)]

            self._cached_all = True

    def add_single_to_cache(self, idx, json):
        """
        Add a single object to cache by json
        """

        if idx not in self.cache:
            simplified = yaml.load(json)
            obj = self.simplifier.build(simplified)

            obj.json = json
            self.index[obj] = idx
            self.cache[idx] = obj

            if self.has_name:
                name = self.storage.variables[self.prefix + '_name'][idx]
                setattr(obj, '_name', name)
                if name != '':
                    self._update_name_in_cache(obj._name, idx)

    def _save(self, obj, idx):
        self.vars['json'][idx] = obj

    @property
    def last(self):
        """
        Returns the last generated trajectory. Useful to continue a run.

        Returns
        -------
        Trajectoy
            the actual trajectory object
        """
        return self.load(len(self) - 1)

    @property
    def first(self):
        """
        Returns the last stored object. Useful to continue a run.

        Returns
        -------
        Object
            the actual last stored object
        """
        return self.load(0)

    def free(self):
        """
        Return the number of the next free index

        Returns
        -------
        index : int
            the number of the next free index in the storage.
            Used to store a new object.
        """
        count = len(self)
        self._free = set([idx for idx in self._free if idx >= count])
        idx = count
        while idx in self._free:
            idx += 1

        return idx

    def reserve_idx(self, idx):
        """
        Locks an idx as used
        """
        self._free.add(idx)

    def _init(self):
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
        self.storage.createDimension(self.prefix, 0)

        if self.has_name:
            self.init_variable("name", 'str',
                               description='A name',
                               chunksizes=tuple([10240]))

        if self.json:
            self.init_variable("json", 'json',
                               description='A json serialized version of the object',
                               chunksizes=tuple([10240]))

    def _restore(self):
        pass

    # ==============================================================================
    # INITIALISATION UTILITY FUNCTIONS
    # ==============================================================================

    def init_variable(self, name, var_type, dimensions=None, **kwargs):
        """
        Create a new variable in the netCDF storage. This is just a helper
        function to structure the code better.

        Parameters
        ==========
        name : str
            The name of the variable to be created
        var_type : str
            The string representing the type of the data stored in the variable.
            Allowed are strings of native python types in which case the variables
            will be treated as python or a string of the form 'numpy.type' which
            will refer to the numpy data types. Numpy is preferred sinec the api
            to netCDF uses numpy and thus it is faster. Possible input strings are
            `int`, `float`, `long`, `str`, `numpy.float32`, `numpy.float64`,
            `numpy.int8`, `numpy.int16`, `numpy.int32`, `numpy.int64`
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
        """

        # add the main dimension to the var_type

        if type(dimensions) is str:
            dimensions = [dimensions]

        if type(dimensions) is int:
            if dimensions == 1:
                dimensions = ['scalar']
            else:
                dimensions = [dimensions]

        if dimensions is None:
            dimensions = (self.prefix,)
        else:
            dimensions = tuple([self.prefix] + list(dimensions))

        self.storage.create_variable(
            self.prefix + '_' + name,
            var_type=var_type,
            dimensions=dimensions,
            **kwargs
        )

    # ==============================================================================
    # COLLECTIVE VARIABLE UTILITY FUNCTIONS
    # ==============================================================================

    @property
    def op_idx(self):
        """
        Returns a function that returns for an object of this storage the idx.
        This can be used to construct order parameters the return the index
        in this storage. Useful for visualization

        Returns
        -------
        function
            the function that reports the index (int) in this store or None if it is not stored
        """

        def idx(obj):
            return self.index.get(obj, None)

        return idx


    # =============================================================================
    # LOAD/SAVE DECORATORS FOR CACHE HANDLING
    # =============================================================================

    def load(self, idx):
        """
        Returns an object from the storage.

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
        """

        if type(idx) is not str and idx < 0:
            return None

        if not hasattr(self, 'cache'):
            return self._load(idx)

        n_idx = idx

        if type(idx) is str:
            # we want to load by name and it was not in cache.
            if self.has_name:
                if idx in self.name_idx:
                    if len(self.name_idx[idx]) > 1:
                        logger.debug('Found name "%s" multiple (%d) times in storage! Loading first!' % (
                            idx, len(self.cache[idx])))

                    n_idx = self.name_idx[idx][0]
                else:
                    # since it is not found in the cache before. Refresh the cache
                    self.update_name_cache()

                    # and give it another shot
                    if idx in self.name_idx:
                        if len(self.name_idx[idx]) > 1:
                            logger.debug('Found name "%s" multiple (%d) times in storage! Loading first!' % (
                                idx, len(self.cache[idx])))

                        n_idx = self.name_idx[idx][0]
                    else:
                        raise ValueError('str "' + idx + '" not found in storage')

        elif type(idx) is not int:
            raise ValueError('indices of type "%s" are not allowed in named storage' % type(idx).__name__)

        # if it is in the cache, return it
        try:
            obj = self.cache[n_idx]
            logger.debug('Found IDX #' + str(idx) + ' in cache. Not loading!')
            return obj

        except KeyError:
            pass

        # turn into python int if it was a numpy int (in some rare cases!)
        n_idx = int(n_idx)

        logger.debug('Calling load object of type ' + self.content_class.__name__ + ' and IDX #' + str(idx))

        if n_idx >= len(self):
            logger.warning('Trying to load from IDX #' + str(n_idx) + ' > number of object ' + str(len(self)))
            return None
        elif n_idx < 0:
            logger.warning('Trying to load negative IDX #' + str(n_idx) + ' < 0. This should never happen!!!')
            raise RuntimeError('Loading of negative int should result in no object. This should never happen!')
        else:
            obj = self._load(idx)

        self.index[obj] = n_idx

        if self.has_name and hasattr(obj, '_name'):
            setattr(obj, '_name',
                    self.storage.variables[self.prefix + '_name'][idx])
            # make sure that you cannot change the name of loaded objects
            obj.fix_name()

        if obj is not None:
            # update cache there might have been a change due to naming
            self.cache[n_idx] = obj

            # finally store the name of a named object in cache
            if self.has_name and obj._name != '':
                self._update_name_in_cache(obj._name, n_idx)

        return obj

    def save(self, obj, idx=None):
        """
        Saves an object to the storage.

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

        if idx is None:
            if obj in self.index:
                # has been saved so quit and do nothing
                return self.index[obj]
            elif type(obj) is LoaderProxy:
                return obj._idx
            else:
                idx = self.free()
        elif type(idx) is str:
            # Not yet supported
            if self.has_name and obj._name_fixed is False:
                obj.name = idx
        else:
            #assume int like
            idx = int(idx)

        self.index[obj] = idx

        # make sure in nested saving that an IDX is not used twice!
        self.reserve_idx(idx)

        logger.debug('Saving ' + str(type(obj)) + ' using IDX #' + str(idx))
        self._save(obj, idx)

        if self.has_name and hasattr(obj, '_name'):
            # logger.debug('Object ' + str(type(obj)) + ' with IDX #' + str(idx))
            # logger.debug(repr(obj))
            # logger.debug("Cleaning up name; currently: " + str(obj._name))
            if obj._name is None:
                # this should not happen!
                logger.debug("Nameable object has not been initialized correctly. Has None in _name")
                raise AttributeError('_name needs to be a string for nameable objects.')

            obj.fix_name()

            self.storage.variables[self.prefix + '_name'][idx] = obj._name

        # store the name in the cache
        if hasattr(self, 'cache'):
            self.cache[idx] = obj
            if self.has_name and obj._name != '':
                # and also the name, if it has one so we can load by
                # name afterwards from cache
                self._update_name_in_cache(obj._name, idx)

        return idx
