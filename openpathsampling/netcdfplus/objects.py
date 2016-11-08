import logging
from uuid import UUID

from cache import MaxCache, Cache, NoCache, WeakLRUCache
from proxy import LoaderProxy
from base import StorableNamedObject, StorableObject

from collections import OrderedDict
from weakref import WeakKeyDictionary, WeakValueDictionary

logger = logging.getLogger(__name__)
init_log = logging.getLogger('openpathsampling.initialization')


class UUIDDict(OrderedDict):
    def __init__(self):
        OrderedDict.__init__(self)

    @staticmethod
    def id(obj):
        if type(obj) is str:
            return UUID(obj)
        elif type(obj) is UUID:
            return obj
        else:
            return obj.__uuid__

    def __getitem__(self, item):
        return OrderedDict.__getitem__(self, self.id(item))

    def __setitem__(self, key, value, **kwargs):
        OrderedDict.__setitem__(self, self.id(key), value)

    def __delitem__(self, key, **kwargs):
        OrderedDict.__delitem__(self, self.id(key))

    def __contains__(self, item):
        return OrderedDict.__contains__(self, self.id(item))

    def get(self, item, default=None):
        return OrderedDict.get(self, self.id(item), default)


class UUIDDictWeak(WeakKeyDictionary):
    def __init__(self):
        WeakKeyDictionary.__init__(self)

    @staticmethod
    def id(obj):
        if type(obj) is str:
            return UUID(obj)
        elif type(obj) is UUID:
            return obj
        else:
            return obj.__uuid__

    def __getitem__(self, item):
        return WeakKeyDictionary.__getitem__(self, self.id(item))

    def __setitem__(self, key, value, **kwargs):
        WeakKeyDictionary.__setitem__(self, self.id(key), value)

    def __delitem__(self, key, **kwargs):
        WeakKeyDictionary.__delitem__(self, self.id(key))

    def __contains__(self, item):
        return WeakKeyDictionary.__contains__(self, self.id(item))

    def get(self, item, default=None):
        return WeakKeyDictionary.get(self, self.id(item), default)


class ObjectStore(StorableNamedObject):
    """
    Base Class for storing complex objects in a netCDF4 file. It holds a
    reference to the store file.`

    Attributes
    ----------
    content_class : :obj:`openpathsampling.netcdfplus.base.StorableObject`
        a reference to the class type to be stored using this Storage. Must be
        subclassed from :obj:`openpathsampling.netcdfplus.base.StorableObject`
    json : string
        if already computed a JSON Serialized string of the object
    cache : :py:class:`openpathsampling.netcdfplus.cache.Cache`
        a dictionary that holds references to all stored elements by index
        or string for named objects. This is only used for cached access
        if caching is not `False`. Must be of type
        :obj:`openpathsampling.netcdfplus.base.StorableObject` or subclassed.

    """

    allowed_types = [
        'int', 'float', 'long', 'str', 'bool',
        'numpy.float32', 'numpy.float64',
        'numpy.int8', 'numpy.inf16', 'numpy.int32', 'numpy.int64',
        'numpy.uint8', 'numpy.uinf16', 'numpy.uint32', 'numpy.uint64',
        'index', 'length', 'uuid'
    ]

    default_store_chunk_size = 250

    class DictDelegator(object):
        def __init__(self, store, dct):
            self.prefix = store.prefix + '_'
            self.dct = dct

        def __getitem__(self, item):
            return self.dct[self.prefix + item]

    def prefix_delegate(self, dct):
        return ObjectStore.DictDelegator(self, dct)

    default_cache = 10000

    def __init__(self, content_class, json=True, nestable=False):
        """

        Parameters
        ----------
        content_class
        json : bool or str `json` or `jsonobj`
            if `False` the store will not create a json variable for
            serialization if `True` the store will use the json pickling to
            store objects and a single storable object will be serialized and
            not referenced. If a string is given the string is taken as the
            variable type of the json variable. Here only two values are
            allowed: `jsonobj` (equivalent to `True`) or `json` which will
            also reference directly given storable objects.

        nestable : bool
            if `True` this marks the content_class to be saved as nested dict
            objects and not a pointing to saved objects. So the saved complex
            object is only stored once and not split into several objects that
            are referenced by each other in a tree-like fashion

        Notes
        -----
        Usually you want caching, but limited. Recommended is to use an LRUCache
        with a reasonable maximum number of objects that depends on the typical
        number of objects to cache and their size

        The class that takes care of storing data in a file is called a
        `Storage`, so the netCDF+ subclassed `Storage` is a storage.
        The classes that know how to load and save an object from the storage
        are called `Store`, like ObjectStore, SampleStore, etc...

        The difference between `json` and `jsonobj` is subtle. Consider
        storing a complex object. Then there are two ways to do that.
        1. `json`: Store a reference to the object (provided) it is stored and
        2. `jsonobj`: serialize the object and only use references for contained
        objects. All inner objects will always be stored using references.
        The only exception is using nestable. Consider objects that contain
        references to objects of the same type, like e.g. operations in an
        equation (2*3 + 3). Each operation represents a value but each
        operation needs values to operate on. To save such an object you have
        again two options:
        1. `nestable=False`. Store all single objects and always reference
        the contained objects. For an equation that would mean to store several
        objects `op1 = plus(op2, 3), op2 = times(2, 3)`. Since this is correct
        though not intuitive you can also use
        2. `nestable=True`. Store all the serialized objects nested into one
        object (string). For our example this corresponds to
        `plus(times(2,3), 3)`.

        """

        super(ObjectStore, self).__init__()
        self._storage = None
        self.content_class = content_class
        self.prefix = None
        self.cache = NoCache()
        self._free = set()
        self._cached_all = False
        self.nestable = nestable
        self._created = False

        # This will not be stored since its information is contained in the
        # dimension names
        self._dimension_prefix_store = None

        self.variables = dict()
        self.vars = dict()
        self.units = dict()

        self.index = None

        self.proxy_index = WeakValueDictionary()

        if json in [True, False, 'json', 'jsonobj']:
            self.json = json
        else:
            raise ValueError(
                'Valid settings for json are only True, False, `json` or '
                '`jsonobj`.')

        if self.content_class is not None \
                and not issubclass(self.content_class, StorableObject):
            raise ValueError(
                'Content class "%s" must be subclassed from StorableObject.' %
                self.content_class.__name__)

        self.fallback_store = None

    def is_created(self):
        return self._created

    def to_dict(self):
        return {
            'content_class': self.content_class,
            'json': self.json,
            'nestable': self.nestable
        }

    def register_fallback(self, store):
        self.fallback_store = store

    def register(self, storage, prefix):
        """
        Associate the object store to a specific storage with a given prefix

        Parameters
        ----------
        storage : :class:`openpathsampling.netcdfplus.NetCDFPlus`
            the storage to be associated with
        prefix : str
            the name under which

        """
        self._storage = storage
        self.prefix = prefix

        self.variables = self.prefix_delegate(self.storage.variables)
        self.units = self.prefix_delegate(self.storage.units)
        self.vars = self.prefix_delegate(self.storage.vars)

        if self.reference_by_uuid:
            self.index = self.create_uuid_index()
        else:
            self.index = self.create_int_index()

    def create_uuid_index(self):
        return UUIDDict()

    def create_int_index(self):
        return UUIDDictWeak()

    @property
    def reference_by_uuid(self):
        return self.storage.reference_by_uuid

    def restore(self):
        if self.reference_by_uuid:
            self.load_indices()

    def load_indices(self):
        uuids = self.vars['uuid'][:]
        for idx, uuid in enumerate(uuids):
            self.index[uuid] = idx

    @property
    def storage(self):
        """Return the associated storage object

        Returns
        -------

        :class:`openpathsampling.netcdfplus.NetCDFPlus`
            the referenced storage object
        """

        if self._storage is None:
            raise RuntimeError(
                'A storage needs to be added to this store to be used! '
                'Use .register() to do so.')

        return self._storage

    def __str__(self):
        return repr(self)

    def __repr__(self):
        return 'store.%s[%s] : %s' % (
            self.prefix,
            self.content_class.__name__ if self.content_class is not None else
            'None/ANY',
            str(len(self)) + ' object(s)' if self._created else
            '(not created)'
        )

    @property
    def simplifier(self):
        """
        Return the simplifier instance used to create JSON serialization

        Returns
        -------
        :class:`openpathsampling.netcdfplus.dictify.StorableObjectJSON`
            the simplifier object used in the associated storage

        """
        return self.storage.simplifier

    def set_caching(self, caching):
        """
        Set the caching mode for this store

        Parameters
        ----------
        caching : :class:`openpathsampling.netcdfplus.Cache`

        """
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
        obj : :class:`openpathsampling.netcdfplus.base.StorableObject`
            the object that can be stored in this store for which its index is
            to be returned

        Returns
        -------
        int or `None`
            The integer index of the given object or `None` if it is not
            stored yet
        """
        return self.index[obj]

    def __iter__(self):
        """
        Add iteration over all elements in the storage
        """
        if self.reference_by_uuid:
            # we want to iterator in the order object were saved!
            for uuid in self.index:
                yield self.load(uuid)
        else:
            for idx in range(len(self)):
                yield self.load(idx)

    def __len__(self):
        """
        Return the number of stored objects

        Returns
        -------
        int
            number of stored objects

        """
        return len(self.storage.dimensions[self.prefix])

    def write(self, variable, idx, obj, attribute=None):
        if attribute is None:
            attribute = variable

        var = self.vars[variable]
        val = getattr(obj, attribute)

        var[int(idx)] = val

        if var.var_type.startswith('lazy'):
            proxy = var.store.proxy(val)
            if isinstance(obj, LoaderProxy):
                # for a loader proxy apply it to the real object
                setattr(obj.__subject__, attribute, proxy)
            else:
                setattr(obj, attribute, proxy)

    def proxy(self, item):
        """
        Return a proxy of a object for this store

        Parameters
        ----------
        item : :py:class:`openpathsampling.netcdfplus.base.StorableObject`
            or int The item or index that points to an object in this store
            and to which a proxy is requested.

        Returns
        -------

        """
        if item is None:
            return None

        tt = type(item)
        if self.reference_by_uuid:
            if tt is int:
                idx = self.vars['uuid'][item]
            elif tt is UUID:
                idx = item
            elif tt in [str, unicode]:
                if item[0] == '-':
                    return None
                idx = UUID(item)
            else:
                idx = item.__uuid__
        else:
            if tt is int:
                if item in self.cache:
                    idx = self.cache[item].__uuid__
                    self.proxy_index[item] = idx
                elif item in self.proxy_index:
                    idx = self.proxy_index[item]
                    self.index[idx] = item
                else:
                    # apparently we want a proxy for a non-existing object
                    # so we create a new UUID and tell the storage
                    # that we associate the UUID with that index
                    idx = StorableObject.get_uuid()
                    self.index[idx] = item
                    self.proxy_index[item] = idx
            else:
                # idx = self.index.get(item)
                idx = item.__uuid__
                if item in self.index:
                    self.proxy_index[self.index[item]] = idx
                else:
                    return item

        return LoaderProxy(self, idx)

    def __getitem__(self, item):
        """
        Enable numpy style selection of object in the store
        """
        try:
            if type(item) is int:
                if item < 0:
                    item += len(self)
                return self.load(item)
            elif type(item) is str or type(item) is UUID:
                return self.load(item)
            elif type(item) is slice:
                return [self.load(idx)
                        for idx in range(*item.indices(len(self)))]
            elif type(item) is list:
                return [self.load(idx) for idx in item]
            elif item is Ellipsis:
                return iter(self)
        except KeyError:
            return None

    def get(self, item):
        try:
            return self[item]
        except KeyError:
            return None

    def _load(self, idx):
        obj = self.vars['json'][idx]
        return obj

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

    def _save(self, obj, idx):
        self.vars['json'][idx] = obj

    @property
    def last(self):
        """
        Returns the last generated trajectory. Useful to continue a run.

        Returns
        -------
        :py:class:`openpathsampling.netcdfplus.base.StorableObject`
            the last stored object in this store
        """
        return self.load(len(self) - 1)

    @property
    def first(self):
        """
        Returns the first stored object.

        Returns
        -------
        :py:class:`openpathsampling.netcdfplus.base.StorableObject`
            the actual first stored object
        """
        return self.load(0)

    def free(self):
        """
        Return the number of the next free index for this store

        Returns
        -------
        index : int
            the number of the next free index in the storage.
            Used to store a new object.
        """

        # start at first free position in the storage
        idx = len(self)

        # and skip also reserved potential stored ones
        while idx in self._free:
            idx += 1

        return idx

    def reserve_idx(self, idx):
        """
        Locks an idx as used

        Parameters
        ----------
        idx : int
            the integer index to be reserved
        """
        self._free.add(idx)

    def release_idx(self, idx):
        """
        Releases a lock on an idx

        Parameters
        ----------
        idx : int
            the integer index to be released
        """
        self._free.discard(idx)

    def initialize(self):
        """
        Initialize the associated storage to allow for object storage. Mainly
        creates an index dimension with the name of the object.
        """
        # define dimensions used for the specific object

        self.storage.create_dimension(self.prefix, 0)

        if self.json:
            jsontype = 'jsonobj'
            if type(self.json) is str:
                jsontype = self.json

            self.create_variable(
                "json",
                jsontype,
                description='A json serialized version of the object',
                chunksizes=tuple([65536])
            )

        if self.storage.reference_by_uuid:
            # TODO: Change to 16byte string
            self.create_variable(
                "uuid", 'uuid',
                description='The uuid of the object',
                chunksizes=tuple([65536])
            )

        self._created = True

    # ==============================================================================
    # INITIALISATION UTILITY FUNCTIONS
    # ==============================================================================

    def create_variable(
            self,
            name,
            var_type,
            dimensions=None,
            chunksizes=None,
            **kwargs
    ):
        """
        Create a new variable in the netCDF storage. This is just a helper
        function to structure the code better.

        Parameters
        ==========
        name : str
            The name of the variable to be created
        var_type : str
            The string representing the type of the data stored in the variable.
            Allowed are strings of native python types in which case the
            variables will be treated as python or a string of the form
            'numpy.type' which will refer to the numpy data types. Numpy is
            preferred since the api to netCDF uses numpy and thus it is faster.
            Possible input strings are `int`, `float`, `long`, `str`,
            `numpy.float32`, `numpy.float64`, `numpy.int8`, `numpy.int16`,
            `numpy.int32`, `numpy.int64`
        dimensions : str or tuple of str
            A tuple representing the dimensions used for the netcdf variable.
            If not specified then the default dimension of the storage is used.
        simtk_units : str
            A string representing the units used if the var_type is `float`
            the units is set to `none`
        description : str
            A string describing the variable in a readable form.
        variable_length : bool
            If true the variable is treated as a variable length (list) of the
            given type. A built-in example for this type is a string which is
            a variable length of char. This make using all the mixed
            stuff superfluous
        chunksizes : tuple of int or int
            A tuple of ints per number of dimensions. This specifies in what
            block sizes a variable is stored. Usually for object related stuff
            we want to store everything of one object at once so this is often
            (1, ..., ...). A single int is interpreted as a tuple with one
            entry.
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

        store_chunk_size = ObjectStore.default_store_chunk_size

        if chunksizes is None and len(dimensions) == 1:
            chunksizes = (store_chunk_size, )
        elif chunksizes is not None and dimensions[-1] == '...' \
                and len(dimensions) == len(chunksizes) + 2:
            chunksizes = tuple([store_chunk_size] + list(chunksizes))
        elif chunksizes is not None and dimensions[-1] != '...' \
                and len(dimensions) == len(chunksizes) + 1:
            chunksizes = tuple([store_chunk_size] + list(chunksizes))

        if self.dimension_prefix:
            dimensions = tuple(
                [dimensions[0]] +
                [
                    self.dimension_prefix + dim if type(dim) is str and
                    dim != '...' else dim for dim in dimensions[1:]
                ]
            )
            chunksizes = tuple(
                [chunksizes[0]] +
                [
                    self.dimension_prefix + chs
                    if type(chs) is str else chs for chs in chunksizes[1:]
                ]
            )

        self.storage.create_variable(
            self.prefix + '_' + name,
            var_type=var_type,
            dimensions=dimensions,
            chunksizes=chunksizes,
            **kwargs
        )

    @property
    def dimension_prefix(self):
        if self._dimension_prefix_store is not None:
            return self._dimension_prefix_store.prefix
        else:
            return ''

    def set_dimension_prefix_store(self, prefix_store=None):
        """
        Select which store or none should be used to prefix dimension names

        If you want to create multiple instances of a store and these should
        have differently long dimensions you need unique names for these. This
        way you can select a store and the dimensions will be prefixed with the
        stores prefix

        Parameters
        ----------
        prefix_store : :obj:`openpathsampling.netcdf.ObjectStore`
            the store from which to use its prefix / name to prefix
            dimension names

        """
        self._dimension_prefix_store = prefix_store

    # ==========================================================================
    # LOAD/SAVE DECORATORS FOR CACHE HANDLING
    # ==========================================================================

    def load(self, idx):
        """
        Returns an object from the storage.

        Parameters
        ----------
        idx : int
            the integer index of the object to be loaded

        Returns
        -------
        :py:class:`openpathsampling.netcdfplus.base.StorableObject`
            the loaded object
        """

        if type(idx) is UUID:
            if idx in self.index:
                n_idx = int(self.index[idx])
            else:
                if self.fallback_store is not None:
                    return self.fallback_store.load(idx)
                elif self.storage.fallback is not None:
                    return self.storage.fallback.stores[self.name].load(idx)
                else:
                    raise ValueError(
                        'str %s not found in storage or fallback' % idx)

        elif type(idx) is not int:
            raise ValueError((
                'indices of type "%s" are not allowed in named storage '
                '(only str and int)') % type(idx).__name__
            )
        else:
            n_idx = int(idx)

        if n_idx < 0:
            return None

        # if it is in the cache, return it
        try:
            obj = self.cache[n_idx]
            logger.debug('Found IDX #' + str(idx) + ' in cache. Not loading!')
            return obj

        except KeyError:
            pass

        logger.debug(
            'Calling load object of type `%s` @ IDX #%d' %
            (self.content_class.__name__, n_idx))

        if n_idx >= len(self):
            logger.warning(
                'Trying to load from IDX #%d > number of object %d' %
                (n_idx, len(self)))
            return None
        elif n_idx < 0:
            logger.warning((
                'Trying to load negative IDX #%d < 0. '
                'This should never happen!!!') % n_idx)
            raise RuntimeError(
                'Loading of negative int should result in no object. '
                'This should never happen!')
        else:
            obj = self._load(n_idx)

        logger.debug(
            'Calling load object of type %s and IDX # %d ... DONE' %
            (self.content_class.__name__, n_idx))

        if obj is not None:
            self._get_id(n_idx, obj)

            # update cache there might have been a change due to naming
            self.index[obj] = n_idx
            self.cache[n_idx] = obj

            logger.debug(
                'Try loading UUID object of type %s and IDX # %d ... DONE' %
                (self.content_class.__name__, n_idx))

        logger.debug(
            'Finished load object of type %s and IDX # %d ... DONE' %
            (self.content_class.__name__, n_idx))

        return obj

    def reference(self, obj):
        if self.reference_by_uuid:
            return obj.__uuid__
        else:
            return self.index.get(obj)

    def remember(self, obj):
        """
        Tell a store that an obj should be assumed as stored

        This is useful, if you do not want to store an object in a specific
        store. Especially to make sure snapshots are not stored multiple times

        Parameters
        ----------
        obj : :py:class:`openpathsampling.netcdfplus.base.StorableObject`
            the object to be fake stored

        """

        if obj not in self.index:
            self.index[obj] = -2

    def forget(self, obj):
        """
        This will revert remembering non-stored objects.

        Stored objects cannot be forgotten

        Parameters
        ----------
        obj : :py:class:`openpathsampling.netcdfplus.base.StorableObject`
            the object to be forgotten

        """

        if obj in self.index:
            if self.index[obj] == -2:
                del self.index[obj]

    def save(self, obj, idx=None):
        """
        Saves an object to the storage.

        Parameters
        ----------
        obj : :class:`openpathsampling.netcdfplus.base.StorableObject`
            the object to be stored
        idx : int or string or `None`
            the index to be used for storing. This is highly discouraged since
            it changes an immutable object (at least in the storage). It is
            better to store also the new object and just ignore the
            previously stored one.

        """

        if obj in self.index:
            # has been saved so quit and do nothing
            if not self.index[obj] == -1:
                return self.reference(obj)

            # numbers other than -1 are reserved for other things

        if isinstance(obj, LoaderProxy):
            if obj._store is self:
                # is a proxy of a saved object so do nothing
                return obj._idx
            else:
                # it is stored but not in this store so we try storing the
                # full snapshot which might be still in cache or memory
                # if that is not the case it will be stored again. This can
                # happen when you load from one store save to another. And load
                # again after some time while the cache has been changed and try
                # to save again the loaded object. We will not explicitly store
                # a table that matches objects between different storages.
                return self.save(obj.__subject__)

        if not isinstance(obj, self.content_class):
            raise ValueError((
                'This store can only store object of base type "%s". Given '
                'obj is of type "%s". You might need to use another store.')
                % (self.content_class, obj.__class__.__name__)
            )

        n_idx = self.free()

        # mark as saved so circular dependencies will not cause infinite loops
        self.index[obj] = n_idx

        # make sure in nested saving that an IDX is not used twice!
        self.reserve_idx(n_idx)

        logger.debug('Saving ' + str(type(obj)) + ' using IDX #' + str(n_idx))

        try:
            self._save(obj, n_idx)

            # store the name in the cache
            if hasattr(self, 'cache'):
                self.cache[n_idx] = obj

        except:
            # in case we did not succeed remove the mark as being saved
            del self.index[obj]
            self.release_idx(n_idx)
            raise

        self.release_idx(n_idx)
        self._set_id(n_idx, obj)

        return self.reference(obj)

    def __setitem__(self, key, value):
        """
        Enable saving using __setitem__

        This only supports writing `store[...] = value`. Not sure if this is
        ever used.

        """
        if key is Ellipsis:
            key = None

        self.save(value, key)

    def load_single(self, idx):
        return self._load(idx)

    def load_range(self, start, end):
        return map(self._load, range(start, end))

    def add_single_to_cache(self, idx, json):
        """
        Add a single object to cache by json

        Parameters
        ----------
        idx : int
            the index where the object was stored
        json : str
            json string the represents a serialized version of the stored object
        """

        if idx not in self.cache:
            obj = self.simplifier.from_json(json)

            self._get_id(idx, obj)

            self.cache[idx] = obj
            self.index[obj] = idx

            return obj

    def uuid(self, uuid):
        """
        Return last object with a given uuid

        Parameters
        ----------
        uuid : str
            the uuid to be searched for

        Returns
        -------
        :py:class:`openpathsampling.netcdfplus.base.StorableObject`
            the last object with a given uuid. This is to mimic an immutable
            object. Once you (re-)save with the same uuid you replace the old
            one and hence you leed to load the last stored one.

        """
        return self.load(uuid)

    def _set_id(self, idx, obj):
        if self.reference_by_uuid:
            self.vars['uuid'][idx] = obj.__uuid__

    def _get_id(self, idx, obj):
        if self.reference_by_uuid:
            obj.__uuid__ = self.vars['uuid'][idx]
        else:
            # check if there exists already a proxy with that idx
            if idx in self.proxy_index:
                obj.__uuid__ = self.proxy_index[idx]


class NamedObjectStore(ObjectStore):
    def __init__(self, content_class, json=True, nestable=False):
        super(NamedObjectStore, self).__init__(
            content_class=content_class,
            json=json,
            nestable=nestable
        )

        self._names_loaded = False
        self._name_idx = dict()

        if self.content_class is not None \
                and not issubclass(self.content_class, StorableNamedObject):
            raise ValueError((
                'Content class "%s" must be subclassed from '
                'StorableNamedObject.') %
                self.content_class.__name__
            )

    def initialize(self):
        """
        Initialize the associated storage to allow for object storage. Mainly
        creates an index dimension with the name of the object.
        """
        super(NamedObjectStore, self).initialize()

        self.create_variable(
            "name", 'str',
            description='The name of the object',
            chunksizes=tuple([65536])
        )

    def add_single_to_cache(self, idx, json):
        """
        Add a single object to cache by json

        Parameters
        ----------
        idx : int
            the index where the object was stored
        json : str
            json string the represents a serialized version of the stored object
        """

        if idx not in self.cache:
            obj = super(NamedObjectStore, self).add_single_to_cache(idx, json)

            name = self.storage.variables[self.prefix + '_name'][idx]
            setattr(obj, '_name', name)
            if name != '':
                self._update_name_in_cache(obj._name, idx)

    @property
    def name_idx(self):
        """
        Returns a dictionary of all names pointing to stored indices

        Returns
        -------
        dict of str : set
            A dictionary that has all stored names as keys and the values are a
            set of indices where an object with this name is found.

        """

        # if not done already cache names once
        if not self._names_loaded:
            self.update_name_cache()

        return self._name_idx

    def update_name_cache(self):
        """
        Update the internal name cache with all stored names in the store.

        This allows to load by name for named objects
        """
        if not self._names_loaded:
            for idx, name in enumerate(
                    self.storage.variables[self.prefix + "_name"][:]):
                self._update_name_in_cache(name, idx)

            self._names_loaded = True

    def _update_name_in_cache(self, name, idx):
        # make sure to cast unicode to str
        name = str(name)
        if name != '':
            if name not in self._name_idx:
                self._name_idx[name] = {idx}
            else:
                if idx not in self._name_idx[name]:
                    self._name_idx[name].add(idx)

    def find(self, name):
        """
        Return last object with a given name

        Parameters
        ----------
        name : str
            the name to be searched for

        Returns
        -------
        :py:class:`openpathsampling.netcdfplus.base.StorableObject`
            the last object with a given name. This is to mimic immutable
            objects. Once you (re-)save with the same name you replace the
            old one and hence you leed to load the last stored one.

        """
        return self.load(name)

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
        return sorted(list(self.name_idx[name]))

    def find_all(self, name):
        if len(self.name_idx[name]) > 0:
            return self[sorted(list(self.name_idx[name]))]

    # ==========================================================================
    # LOAD/SAVE DECORATORS FOR CACHE HANDLING
    # ==========================================================================

    def load(self, idx):
        """
        Returns an object from the storage.

        Parameters
        ----------
        idx : int or str
            either the integer index of the object to be loaded or a string
            (name) for named objects. This will always return the last object
            found with the specified name. This allows to effectively change
            existing objects.

        Returns
        -------
        :py:class:`openpathsampling.netcdfplus.base.StorableNamedObject`
            the loaded object
        """

        if type(idx) is int and idx < 0:
            return None

        n_idx = idx

        if type(idx) is str:
            # we want to load by name and it was not in cache.
            if idx in self.name_idx:
                if len(self.name_idx[idx]) > 1:
                    logger.debug((
                        'Found name "%s" multiple (%d) times in storage! '
                        'Loading last!') % (
                        idx, len(self.cache[idx])))

                n_idx = sorted(list(self.name_idx[idx]))[-1]
            else:
                raise ValueError('str "' + idx + '" not found in storage')

        elif type(idx) is UUID:
            pass

        elif type(idx) is not int:
            raise ValueError((
                'indices of type "%s" are not allowed in named '
                'storage (only str and int)') %
                type(idx).__name__
            )

        obj = super(NamedObjectStore, self).load(n_idx)

        if obj is not None:
            n_idx = self.index[obj]
            setattr(obj, '_name',
                    self.storage.variables[self.prefix + '_name'][n_idx])
            # make sure that you cannot change the name of loaded objects
            obj.fix_name()

            # finally store the name of a named object in cache
            self._update_name_in_cache(obj._name, n_idx)

        return obj

    def save(self, obj, idx=None):
        """
        Saves an object to the storage.

        Parameters
        ----------
        obj : :py:class:`openpathsampling.netcdfplus.base.StorableNamedObject`
            the object to be stored
        idx : int or string or `None`
            the index to be used for storing. This is highly discouraged since
            it changes an immutable object (at least in the storage). It is
            better to store also the new object and just ignore the
            previously stored one.

        """

        is_str = type(idx) is str

        if not is_str and idx is not None:
            raise ValueError(
                'Unsupported index type (only str or None allowed).')

        name = obj._name

        if is_str:
            obj.name = idx
            name = obj._name

        if name is None:
            # this should not happen!
            logger.debug(
                "Nameable object has not been initialized correctly. "
                "Has None in _name")
            raise AttributeError(
                '_name needs to be a string for nameable objects.')

        obj_fixed = obj._name_fixed
        obj_name = obj._name
        # we fix the name just in case we try in recursive saving to store
        # one object twice with different names. Storing with the same name
        # is fine!
        obj.fix_name()

        try:
            reference = super(NamedObjectStore, self).save(obj)
        except:
            # if saving did not work unlock the name if is was un-fixed before
            obj._name_fixed = obj_fixed
            obj._name = obj_name
            raise

        n_idx = self.index[obj]
        self.storage.variables[self.prefix + '_name'][n_idx] = name
        self._update_name_in_cache(name, n_idx)

        return reference


class UniqueNamedObjectStore(NamedObjectStore):

    # ==========================================================================
    # LOAD/SAVE DECORATORS FOR CACHE HANDLING
    # ==========================================================================

    def __init__(self, content_class, json=True, nestable=False):
        super(UniqueNamedObjectStore, self).__init__(
            content_class=content_class,
            json=json,
            nestable=nestable)

        self._free_name = set()

    def reserve_name(self, name):
        """
        Locks a name as used

        Parameters
        ----------
        name : str
            the name to be locked for storage
        """
        if name != "":
            self._free_name.add(name)

    def release_name(self, name):
        """
        Releases a locked name

        Parameters
        ----------
        name : str
            the name to be released for being used as a name
        """
        self._free_name.discard(name)

    def is_name_locked(self, name):
        """
        Test whether in a unique name store a name is already taken

        Parameters
        ----------
        name : str or `None`
            the name to be tested.

        Returns
        -------
        bool
            the result of the test. If the name exists or is reserved during
            a saving event this will return `True` and return `False` if the
            name is free.

        """
        if name is None:
            return False

        return name in self.name_idx or name in self._free_name

    def save(self, obj, idx=None):
        """
        Saves an object to the storage.

        Parameters
        ----------
        obj : :py:class:`openpathsampling.netcdfplus.base.StorableNamedObject`
            the object to be stored
        idx : string or `None`
            the index to be used for storing. This is highly discouraged since
            it changes an immutable object (at least in the storage). It is
            better to store also the new object and just ignore the
            previously stored one.

        """

        is_str = type(idx) is str

        if not is_str and idx is not None:
            raise ValueError(
                'Unsupported index type (only str or None allowed).')

        name = obj._name
        fixed = obj._name_fixed
        err = list()

        if is_str:
            if fixed:
                if name != idx:
                    # saving fixed under different name is not possible.
                    if obj in self:
                        err.append(
                            ('Cannot rename object to "%s". '
                             'Already saved with name "%s" !') % (idx, name)
                        )
                    else:
                        err.append(
                            ('Cannot rename object to "%s". '
                             'Already fixed name "%s" !') % (idx, name)
                        )

                        if self.is_name_locked(name):
                            err.append(
                                ('Current name "%s" is also already taken in '
                                 'unique name store. This means you cannot '
                                 'save object "%s" at all. '
                                 'In general this should not happen to unsaved '
                                 'objects unless you fixed the name of the '
                                 'object yourself. Check your code '
                                 'for the generation of objects of the same '
                                 'name.') %
                                (name, obj)
                            )
                        else:
                            err.append(
                                ('Current name "%s" is still free. Saving '
                                 'without giving a specific name '
                                 'should work. If that is what you want '
                                 'to do.') % name
                            )
                else:
                    # already fixed, but with same name. Okay. Check if stored
                    if obj in self.index:
                        return self.reference(obj)
            else:
                # name is not fixed yet. Check, if we can save or name is taken
                if self.is_name_locked(idx):
                    err.append(
                        ('New name "%s" already taken in unique name store. ' +
                         'Try different name instead.') % idx
                    )

                    if self.is_name_locked(name):
                        err.append((
                            'Current name "%s" already taken in unique '
                            'name store. ') % name
                        )
                    else:
                        err.append(
                            ('Current name "%s" is still free. Saving without '
                             'giving a specific name should work') % name
                        )
        else:
            if fixed:
                # no new name, but fixed. Check if already stored.
                if obj in self.index:
                    return self.reference(obj)

                # if not stored yet check if we could
                if self.is_name_locked(name):
                    err.append(
                        ('Current name "%s" is already taken in unique name '
                         'store. This means you cannot save object "%s" at '
                         'all. In general this should not happen to unsaved '
                         'objects unless you fixed the name of the object '
                         'yourself. Check your code for the generation of '
                         'objects of the same name.') %
                        (name, obj)
                    )
            else:
                # no new name and not fixed. Just check if current name is taken
                if self.is_name_locked(name):
                    err.append(
                        ('Current name "%s" is already taken in unique name '
                         'store %s. Try renaming object or saving using other '
                         'name.') % (name, self.name)
                    )

        if len(err) > 0:
            raise RuntimeWarning('/n'.join(err))

        # no errors, reserve the name for nested saving and actually call save
        self.reserve_name(name)

        try:
            reference = super(UniqueNamedObjectStore, self).save(obj, idx)
        finally:
            self.release_name(name)

        return reference


class VariableStore(ObjectStore):
    def __init__(self, content_class, var_names):
        super(VariableStore, self).__init__(
            content_class,
            json=False
        )

        # TODO: determine var_names automatically from content_class
        # problem is that some decorators, e.g. using delayed loader
        # hide the actual __init__ signature and so we cannot determine
        # what variables to store. Could be 2.0

        if not issubclass(content_class, StorableObject):
            raise ValueError(('Content_class %s must be subclassed from '
                             'StorableObject') % content_class.__name__)

        self.var_names = var_names
        self._cached_all = False

    def to_dict(self):
        return {
            'content_class': self.content_class,
            'var_names': self.var_names
        }

    def _save(self, obj, idx):
        for var in self.var_names:
            self.write(var, idx, obj)

    def _load(self, idx):
        # attr = {var: self.vars[var][idx] for var in self.var_names}
        args = [ self.vars[var][idx] for var in self.var_names]
        return self.content_class(*args)

    def initialize(self):
        super(VariableStore, self).initialize()

        # Add here the stores to be imported
        # self.create_variable('name', 'var_type')

    def all(self):
        self.cache_all()
        return self

    def cache_all(self, part=None):
        """Load all samples as fast as possible into the cache

        Parameters
        ----------
        part : list of int or `None`
            If `None` (default) all samples will be loaded. Otherwise the
            list of indices in `part` will be loaded into the cache

        """
        max_length = self.cache.size[0]
        max_length = len(self) if max_length < 0 else max_length

        if part is None:
            length = min(len(self), max_length)
            part = range(length)
        else:
            part = sorted(list(set(part())))
            length = min(len(part), max_length)
            part = part[:length]

        if not part:
            return

        # just in case we saved the var_names in another order and so we are
        # backwards compatible
        var_names = self.content_class.args()[1:]

        if not self._cached_all:
            data = zip(*[
                self.vars[var][part]
                for var in var_names
            ])

            [self.add_to_cache(idx, v) for idx, v in zip(part, data)]

            self._cached_all = True

    def add_to_cache(self, idx, data):
        if idx not in self.cache:
            # attr = {var: self.vars[var].getter(data[nn])
            #         for nn, var in enumerate(self.var_names)}
            obj = self.content_class(*data)
            self._get_id(idx, obj)

            self.index[obj] = idx
            self.cache[idx] = obj


class DictStore(NamedObjectStore):
    def __init__(self):
        super(DictStore, self).__init__(
            None,
            json='json'
        )

    def to_dict(self):
        return {}

    def load(self, idx):
        """
        Returns an object from the storage.

        Parameters
        ----------
        idx : str
            a string (name) of the objects. This will always return the
            last object found with the specified name. If immutable is true
            for the store it assures that there is only a single object per name

        Returns
        -------
        :class:`openpathsampling.netcdfplus.base.StorableObject`
            the loaded object
        """

        if type(idx) is str:
            n_idx = -1

            # we want to load by name and it was not in cache.
            if idx not in self.name_idx:
                logger.debug('Name "%s" not found in the storage!' % idx)
                raise KeyError('str "' + idx + '" not found in storage')

            if idx in self.name_idx:
                if len(self.name_idx[idx]) > 1:
                    logger.debug((
                        'Found name "%s" multiple (%d) times in storage! '
                        'Loading last!') % (
                        idx, len(self.name_idx[idx])))

                n_idx = sorted(list(self.name_idx[idx]))[-1]

                logger.debug(
                    'Found name "%s" in storage! Loading IDX %d!' %
                    (idx, n_idx))

        elif type(idx) is int:
            n_idx = idx
        else:
            raise ValueError(
                'Unsupported index type (only str and int allowed).')

        # turn into python int if it was a numpy int (in some rare cases!)
        n_idx = int(n_idx)

        logger.debug(
            'Calling load object of type `%s` and @ IDX #%d' %
            (str(self.content_class), n_idx))

        if n_idx >= len(self):
            logger.warning(
                'Trying to load from IDX #' + str(n_idx) +
                ' > number of objects ' + str(len(self))
            )
            raise RuntimeError(
                'Loading of too large int should be attempted. '
                'Problem in name cache. This should never happen!')
        elif n_idx < 0:
            logger.warning(
                'Trying to load negative IDX #' + str(n_idx) + ' < 0. '
                'This should never happen!!!'
            )
            raise RuntimeError(
                'Loading of negative int should result in no object. '
                'This should never happen!'
            )
        else:

            logger.debug('Loading named object from index IDX # %d' % n_idx)
            obj = self._load(n_idx)

            logger.debug(
                'Loading named object from index IDX # %d.. DONE' % n_idx)

        return obj

    def restore(self):
        pass

    def save(self, obj, idx=None):
        """
        Saves an object to the storage.

        Parameters
        ----------
        obj : :py:class:`openpathsampling.netcdfplus.base.StorableObject`
            the object to be stored
        idx : string or `None`
            the string index to be used for storing. Objects will not be
            replaced but stored again with the same name. When loading the
            last stored object under the idx is retrieved. Effectively mimicking
            a mutual dict with versioning. We usually encourage for most cases
            to use the immutual dict class
            :class:`openpathsampling.netcdfplus.ImmutableDictStore` instead
            to avoid ambiguity in stored objects.

        See Also
        --------
        :class:`openpathsampling.netcdfplus.ImmutableDictStore`

        """

        if idx is None:
            # a DictStore needs a specific name
            raise ValueError(
                'Saving in a DictStore without specifying a '
                'string key is not allowed. ')

        if type(idx) is not str:
            # key needs to be a string
            raise ValueError(
                'Index `%s` for DictStore needs to be a string! ' % idx)

        n_idx = int(self.free())
        # make sure in nested saving that an IDX is not used twice!
        self.reserve_idx(n_idx)

        logger.debug(
            'Saving `%s` with name `%s` @ IDX #%d' %
            (str(obj.__class__), idx, n_idx))
        self._save(obj, n_idx)

        self.storage.variables[self.prefix + '_name'][n_idx] = idx
        self._update_name_in_cache(idx, n_idx)

        return n_idx

    def keys(self):
        return self.name_idx.keys()

    def iterkeys(self):
        return self.name_idx.iterkeys()

    def __iter__(self):
        return self.iterkeys()

    def iteritems(self):
        for name in self:
            yield name, self[name]

    def get(self, idx, default=None):
        try:
            return self.load(idx)
        except KeyError:
            return default


class ImmutableDictStore(DictStore):

    def save(self, obj, idx=None):
        """
        Saves an object to the storage.

        Parameters
        ----------
        obj : :py:class:`openpathsampling.netcdfplus.base.StorableObject`
            the object to be stored
        idx : int or string or `None`
            the index to be used for storing. This is highly discouraged since
            it changes an immutable object (at least in the storage). It is
            better to store also the new object and just ignore the
            previously stored one.

        """

        if idx in self.name_idx:
            # immutable means no duplicates, so quit
            raise RuntimeWarning(
                'Cannot re-save existing key "%s" in '
                'immutable dict store.' % idx
            )

        return super(ImmutableDictStore, self).save(obj, idx)


class IndexedObjectStore(ObjectStore):
    """
    ObjectStore storing objects in arbitrary order

    This has a prefilled .index which knows at which position a certain
    index is stored. This way you can circumvent holes and keep the file smaller
    """

    # ==========================================================================
    # LOAD/SAVE DECORATORS FOR CACHE HANDLING
    # ==========================================================================

    def load(self, idx):
        """
        Returns an object from the storage.

        Parameters
        ----------
        idx : int
            the integer index of the object to be loaded

        Returns
        -------
        :py:class:`openpathsampling.netcdfplus.base.StorableObject`
            the loaded object
        """

        # we want to load by uuid and it was not in cache.
        if idx in self.index:
            n_idx = self.index[idx]
        else:
            raise KeyError(idx)

        if n_idx < 0:
            return None

        # if it is in the cache, return it
        try:
            obj = self.cache[n_idx]
            return obj

        except KeyError:
            pass

        obj = self._load(n_idx)
        self.cache[n_idx] = obj

        return obj

    @property
    def reference_by_uuid(self):
        return False

    def create_int_index(self):
        return dict()

    def save(self, obj, idx=None):
        """
        Saves an object to the storage.

        Parameters
        ----------
        obj : :py:class:`openpathsampling.netcdfplus.base.StorableObject`
            the object to be stored
        idx : int or string or `None`
            the index to be used for storing. This is highly discouraged since
            it changes an immutable object (at least in the storage). It is
            better to store also the new object and just ignore the
            previously stored one.

        """

        if idx in self.index:
            # has been saved so quit and do nothing
            return idx

        n_idx = self.free()

        # mark as saved so circular dependencies will not cause infinite loops
        self.index[idx] = n_idx

        # make sure in nested saving that an IDX is not used twice!
        self.reserve_idx(n_idx)

        logger.debug('Saving ' + str(type(obj)) + ' using IDX #' + str(n_idx))

        try:
            self._save(obj, n_idx)
            self.vars['index'][n_idx] = idx

            # store the name in the cache
            if hasattr(self, 'cache'):
                self.cache[n_idx] = obj

        except:
            logger.debug('Problem saving %d !' % n_idx)
            # in case we did not succeed remove the mark as being saved
            del self.index[idx]
            self.release_idx(n_idx)
            raise

        self.release_idx(n_idx)
        self._set_id(n_idx, obj)

        return idx

    def restore(self):
        for pos, idx in enumerate(self.vars['index'][:]):
            self.index[idx] = pos

    def initialize(self):
        super(IndexedObjectStore, self).initialize()

        self.create_variable(
            'index',
            'index'
        )
