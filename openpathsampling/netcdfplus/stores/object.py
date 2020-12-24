import logging
# from uuid import UUID
from weakref import WeakValueDictionary

from openpathsampling.netcdfplus.base import StorableNamedObject, StorableObject
from openpathsampling.netcdfplus.cache import MaxCache, Cache, NoCache, \
    WeakLRUCache
from openpathsampling.netcdfplus.proxy import LoaderProxy

from future.utils import iteritems

import sys
if sys.version_info > (3, ):
    long = int
    unicode = str

logger = logging.getLogger(__name__)
init_log = logging.getLogger('openpathsampling.initialization')


class HashedList(dict):
    def __init__(self):
        super(HashedList, self).__init__()
        dict.__init__(self)
        self._list = []

    def append(self, key):
        dict.__setitem__(self, key, len(self))
        self._list.append(key)

    # noinspection PyCallByClass
    def extend(self, t):
        l = len(self)
        dict.update(self, zip(t, range(l, l + len(t))))
        self._list.extend(t)

    def __setitem__(self, key, value):
        dict.__setitem__(self, key, value)
        self._list[value] = key

    def __getitem__(self, key):
        return dict.__getitem__(self, key)

    def index(self, key):
        return self._list[key]

    def mark(self, key):
        if key not in self:
            dict.__setitem__(self, key, -2)

    def unmark(self, key):
        if key in self:
            dict.__delitem__(self, key)

    def clear(self):
        dict.clear(self)
        self._list = []

    @property
    def list(self):
        return self._list


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
    _restore_non_initial_attr = False

    allowed_types = [
        'int', 'float', 'long', 'str', 'bool',
        'numpy.float32', 'numpy.float64',
        'numpy.int8', 'numpy.inf16', 'numpy.int32', 'numpy.int64',
        'numpy.uint8', 'numpy.uinf16', 'numpy.uint32', 'numpy.uint64',
        'index', 'length', 'uuid'
    ]

    default_store_chunk_size = 256

    _log_debug = False

    class DictDelegator(object):
        def __init__(self, store, dct):
            self.prefix = store.prefix + '_'
            self.dct = dct

        def __getitem__(self, item):
            return self.dct[self.prefix + item]

        def __contains__(self, item):
            return (self.prefix + item) in self.dct

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

        self.attribute_list = {}
        self.cv = {}

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

        self.index = self.create_uuid_index()

    def create_uuid_index(self):
        return HashedList()

    def restore(self):
        self.load_indices()

    def load_indices(self):
        self.index.clear()
        self.index.extend(self.vars['uuid'][:])

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
        return self.index[obj.__uuid__]

    def __iter__(self):
        """
        Add iteration over all elements in the storage
        """
        # we want to iterator in the order object were saved!
        for uuid in self.index._list:
            yield self.load(uuid)

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

        try:
            idx = item.__uuid__
        except AttributeError:
            idx = item

        # tt = type(item)
        # if tt is int:
        #     idx = self.vars['uuid'][item]
        # elif tt is long:
        #     idx = item
        # elif tt in [str, unicode]:
        #     if item[0] == '-':
        #         return None
        #     idx = int(UUID(item))
        # else:
        #
        return LoaderProxy.new(self, idx)

    def __contains__(self, item):
        if item.__uuid__ in self.index:
            return True

        if self.fallback_store is not None and item in self.fallback_store:
            return True

        if self.storage.fallback is not None and item in self.storage.fallback:
            return True

        return False

    def __getitem__(self, item):
        """
        Enable numpy style selection of object in the store
        """
        try:
            if isinstance(item, (long, int)):
                if item < 0:
                    item += len(self)
                return self.load(item)
            elif type(item) is str:
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
        """Clear the cache and force reloading"""

        self.cache.clear()
        self._cached_all = False

    def cache_all(self):
        """Load all samples as fast as possible into the cache"""
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

        # # and skip also reserved potential stored ones
        # while idx in self._free:
        #     idx += 1

        return idx

    # def reserve_idx(self, idx):
    #     """
    #     Locks an idx as used
    #
    #     Parameters
    #     ----------
    #     idx : int
    #         the integer index to be reserved
    #     """
    #     self._free.add(idx)
    #
    # def release_idx(self, idx):
    #     """
    #     Releases a lock on an idx
    #
    #     Parameters
    #     ----------
    #     idx : int
    #         the integer index to be released
    #     """
    #     self._free.discard(idx)

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

        # TODO: Change to 16byte string
        self.create_variable(
            "uuid", 'uuid',
            description='The uuid of the object',
            chunksizes=tuple([65536])
        )

        self._created = True

    # ==========================================================================
    # INITIALISATION UTILITY FUNCTIONS
    # ==========================================================================

    def create_variable(
            self,
            var_name,
            var_type,
            dimensions=None,
            chunksizes=None,
            description=None,
            simtk_unit=None,
            maskable=False
    ):
        """
        Create a new variable in the netCDF storage. This is just a helper
        function to structure the code better.

        Parameters
        ==========
        var_name : str
            The var_name of the variable to be created
        var_type : str
            The string representing the type of the data stored in the
            variable.  Allowed are strings of native python types in which
            case the variables will be treated as python or a string of the
            form 'numpy.type' which will refer to the numpy data types.
            Numpy is preferred sinec the api to netCDF uses numpy and thus
            it is faster. Possible input strings are
            `int`, `float`, `long`, `str`, `numpy.float32`, `numpy.float64`,
            `numpy.int8`, `numpy.int16`, `numpy.int32`, `numpy.int64`, `json`,
            `obj.<store>`, `lazyobj.<store>`
        dimensions : str or tuple of str
            A tuple representing the dimensions used for the netcdf variable.
            If not specified then the default dimension of the storage is used.
            If the last dimension is `'...'` then it is assumed that the
            objects are of variable length. In netCDF this is usually
            referred to as a VLType.  We will treat is just as another
            dimension, but it can only be the last dimension.
        description : str
            A string describing the variable in a readable form.
        chunksizes : tuple of int
            A tuple of ints per number of dimensions. This specifies in what
            block sizes a variable is stored. Usually for object related stuff
            we want to store everything of one object at once so this is often
            (1, ..., ...)
        simtk_unit : str
            A string representing the units used for this variable. Can be
            used with all var_types although it makes sense only for numeric
            ones.
        maskable : bool, default: False
            If set to `True` the values in this variable can only partially
            exist and if they have not yet been written they are filled with
            a fill_value which is treated as a non-set variable. The created
            variable will interpret this values as `None` when returned
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
            self.prefix + '_' + var_name,
            var_type=var_type,
            dimensions=dimensions,
            chunksizes=chunksizes,
            description=description,
            simtk_unit=simtk_unit,
            maskable=maskable
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

        if isinstance(idx, (long, int)):
            if idx < 1000000000:
                n_idx = idx
            elif idx in self.index:
                n_idx = self.index[idx]
            else:
                if self.fallback_store is not None:
                    return self.fallback_store.load(idx)
                elif self.storage.fallback is not None:
                    return self.storage.fallback.stores[self.name].load(idx)
                else:
                    raise ValueError(
                        'str %s not found in storage or fallback' % idx)

        else:
            raise ValueError(
                'indices need to be a 32-byte UUID in long format or a simple int ')

        if n_idx < 0:
            return None

        # if it is in the cache, return it
        try:
            obj = self.cache[n_idx]
            if self._log_debug:
                logger.debug(
                    'Found IDX #' + str(idx) + ' in cache. Not loading!')
            return obj

        except KeyError:
            pass

        if self._log_debug:
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

        if self._log_debug:
            logger.debug(
                'Calling load object of type %s and IDX # %d ... DONE' %
                (self.content_class.__name__, n_idx))

        if obj is not None:
            self._get_id(n_idx, obj)

            # update cache there might have been a change due to naming
            self.cache[n_idx] = obj

            if self._log_debug:
                logger.debug(
                    'Try loading UUID object of type %s and IDX # %d ... DONE' %
                    (self.content_class.__name__, n_idx))

        if self._log_debug:
            logger.debug(
                'Finished load object of type %s and IDX # %d ... DONE' %
                (self.content_class.__name__, n_idx))

        return obj

    @staticmethod
    def reference(obj):
        return obj.__uuid__

    def remember(self, obj):
        """
        Tell a store that an obj should be assumed as stored

        This is useful, if you do not want to store an object in a specific
        store. Especially to make sure attributes are not stored multiple times

        Parameters
        ----------
        obj : :py:class:`openpathsampling.netcdfplus.base.StorableObject`
            the object to be fake stored

        """
        self.index.mark(obj.__uuid__)

    def forget(self, obj):
        """
        This will revert remembering non-stored objects.

        Stored objects cannot be forgotten

        Parameters
        ----------
        obj : :py:class:`openpathsampling.netcdfplus.base.StorableObject`
            the object to be forgotten

        """

        self.index.unmark(obj.__uuid__)

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
        uuid = obj.__uuid__

        if uuid in self.index:
            # has been saved so quit and do nothing
            if not self.index[uuid] == -1:
                return self.reference(obj)

            # numbers other than -1 are reserved for other things

        if isinstance(obj, LoaderProxy):
            if obj._store is self:
                # is a proxy of a saved object so do nothing
                return uuid
            else:
                # it is stored but not in this store so we try storing the
                # full attribute which might be still in cache or memory
                # if that is not the case it will be stored again. This can
                # happen when you load from one store save to another. And load
                # again after some time while the cache has been changed and try
                # to save again the loaded object. We will not explicitly store
                # a table that matches objects between different storages.
                return self.save(obj.__subject__)

        if self.fallback_store is not None and \
                self.storage.exclude_from_fallback:
            if obj in self.fallback_store:
                return self.reference(obj)

        elif self.storage.fallback is not None and \
                self.storage.exclude_from_fallback:
            if obj in self.storage.fallback:
                return self.reference(obj)

        if not isinstance(obj, self.content_class):
            raise ValueError((
                'This store can only store object of base type "%s". Given '
                'obj is of type "%s". You might need to use another store.')
                % (self.content_class, obj.__class__.__name__)
            )

        # n_idx = self.free()
        n_idx = len(self.index)

        # mark as saved so circular dependencies will not cause infinite loops
        self.index.append(uuid)

        # make sure in nested saving that an IDX is not used twice!
        # self.reserve_idx(n_idx)

        logger.debug('Saving ' + str(type(obj)) + ' using IDX #' + str(n_idx))

        try:
            self._save(obj, n_idx)
            self._auto_complete(obj, n_idx)
            self.cache[n_idx] = obj

        except:
            # in case we did not succeed remove the mark as being saved
            del self.index[uuid]
            raise

        # self.release_idx(n_idx)
        self._set_id(n_idx, obj)

        return self.reference(obj)

    def __setitem__(self, key, value):
        """
        Enable saving using __setitem__

        """
        self.save(value, key)

    # def load_single(self, idx):
    #     return self._load(idx)
    #
    # def load_range(self, start, end):
    #     return map(self._load, range(start, end))

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
            self.index[obj.__uuid__] = idx

            return obj

    # def uuid(self, uuid):
    #     """
    #     Return last object with a given uuid
    #
    #     Parameters
    #     ----------
    #     uuid : str
    #         the uuid to be searched for
    #
    #     Returns
    #     -------
    #     :py:class:`openpathsampling.netcdfplus.base.StorableObject`
    #         the last object with a given uuid. This is to mimic an immutable
    #         object. Once you (re-)save with the same uuid you replace the old
    #         one and hence you leed to load the last stored one.
    #
    #     """
    #     return self.load(uuid)

    def _set_id(self, idx, obj):
        self.vars['uuid'][idx] = obj.__uuid__

    def _get_id(self, idx, obj):
        obj.__uuid__ = self.index.index(int(idx))

    # CV SUPPORT

    def _auto_complete(self, obj, pos):
        for attribute, attribute_store in self.attribute_list.items():
            if not attribute_store.allow_incomplete:
                # value = attribute._cache_dict._get(obj)
                # if value is None:
                #     # not in cache so compute it if possible
                #     if attribute._eval_dict:
                #         value = attribute._eval_dict([obj])[0]

                value = attribute(obj)

                if value is not None:
                    if attribute_store.allow_incomplete:
                        attribute_store[obj] = value
                    else:
                        n_idx = pos
                        attribute_store.vars['value'][n_idx] = value
                        attribute_store.cache[n_idx] = value

    def complete_attribute(self, attribute):
        """
        Compute all missing values of a CV and store them


        Parameters
        ----------
        attribute : :obj:`openpathsampling.netcdfplus.PseudoAttribute`


        """
        if attribute not in self.attribute_list:
            return

        attribute_store = self.attribute_list[attribute]
        key_store = self.storage.attributes.key_store(attribute)

        if attribute_store.allow_incomplete:
            # for complete this does not make sense

            # TODO: Make better looping over this to not have
            # to load all the indices at once
            # can be problematic for 10M+ stored attributes
            indices = self.vars['uuid'][:]

            for pos, idx in enumerate(indices):
                if pos not in attribute_store.index:
                    # this value is not stored to go ahead

                    proxy = LoaderProxy.new(key_store, idx)

                    # # get from cache first, this is fastest
                    # value = attribute._cache_dict._get(proxy)
                    #
                    # if value is None:
                    #     # not in cache so compute it if possible
                    #     if attribute._eval_dict:
                    #         value = attribute._eval_dict([proxy])[0]
                    #     else:
                    #         value = None

                    value = attribute(proxy)

                    if value is not None:
                        n_idx = attribute_store.free()

                        attribute_store.vars['value'][n_idx] = value
                        attribute_store.vars['index'][n_idx] = pos
                        attribute_store.index[pos] = n_idx
                        attribute_store.cache[n_idx] = value

    def sync_attribute(self, attribute):
        """
        Store all cached values of a CV in the diskcache

        Parameters
        ----------
        attribute : :obj:`openpathsampling.CollectiveVariable`


        """

        if attribute not in self.attribute_list:
            return

        attribute_store = self.attribute_list[attribute]

        # for complete this does not make sense
        if attribute_store.allow_incomplete:

            # loop all objects in the fast CV cache
            for obj, value in iteritems(attribute._cache_dict.cache):
                if value is not None:
                    pos = self.pos(obj)

                    # if the attribute is not saved, there is nothing we can do
                    if pos is None:
                        continue

                    # if the value is stored, skip it
                    if pos in attribute_store.index:
                        continue

                    n_idx = attribute_store.free()

                    attribute_store.vars['value'][n_idx] = value
                    attribute_store.vars['index'][n_idx] = pos
                    attribute_store.index[pos] = n_idx
                    attribute_store.cache[n_idx] = value

    @staticmethod
    def _get_attribute_name(attribute_idx):
        return 'attribute' + str(attribute_idx)

    def pos(self, obj):
        return self.index.get(obj.__uuid__)

    def pos_uuid(self, uid):
        return self.index.get(uid)

    def add_attribute(
            self, store_cls, attribute, template,
            allow_incomplete=None, chunksize=None):
        """

        Parameters
        ----------
        store_cls : :obj:`openpathsampling.netcdfplus.ValueStore`
        attribute : :obj:`openpathsampling.CollectiveVariable`
        template : :obj:`openpathsampling.engines.Baseattribute`
        chunksize : int
        allow_incomplete : bool

        Returns
        -------
        :obj:`openpathsampling.netcdfplus.ObjectStore`
        int
        """
        if attribute in self.attribute_list:
            return self.attribute_list[attribute]

        key_store = self.storage.attributes.key_store(attribute)

        if allow_incomplete is None:
            allow_incomplete = attribute.diskcache_allow_incomplete
        if chunksize is None:
            chunksize = attribute.diskcache_chunksize

        if template is None:
            template = attribute.diskcache_template

        if not allow_incomplete:
            # in complete mode we force chunk size one to match it to attributes
            # chunksize = self.default_store_chunk_size
            chunksize = self.variables['uuid'].chunking()[0]

        # determine value type and shape
        params = self.storage.get_value_parameters(attribute(template))
        shape = params['dimensions']

        if shape is None:
            chunksizes = None
        else:
            chunksizes = tuple(params['dimensions'])

        # attribute_idx = self.storage.attributes.index[attribute.__uuid__]
        value_store = store_cls(
            attribute.key_class,
            allow_incomplete=allow_incomplete,
            chunksize=chunksize
        )

        store_name = self.name + '_' + attribute.name

        self.storage.create_store(store_name, value_store, False)

        if value_store.allow_incomplete:
            # we are not using the .initialize function here since we
            # only have one variable and only here know its shape
            self.storage.create_dimension(value_store.prefix, 0)

            if shape is not None:
                shape = tuple(list(shape))
                chunksizes = tuple([chunksize] + list(chunksizes))
            else:
                shape = tuple()
                chunksizes = tuple([chunksize])

            # create the variable
            value_store.create_variable(
                'value',
                var_type=params['var_type'],
                dimensions=shape,
                chunksizes=chunksizes,
                simtk_unit=params['simtk_unit'],
            )

            value_store.create_variable('index', 'index')

        else:
            # todo: seems to be a bug in NetCDF4. Need to set chunksize to 1
            # see Issue https://github.com/Unidata/netcdf4-python/issues/566
            # I assume this will still work as expected.

            # chunksize = self.default_store_chunk_size
            # chunksize = self.variables['uuid'].chunking()[0]
            chunksize = 1
            if shape is not None:
                shape = tuple([self.name] + list(shape))
                chunksizes = tuple([chunksize] + list(chunksizes))
            else:
                shape = tuple([self.name])
                chunksizes = tuple([chunksize])

            # create the variable
            value_store.storage.create_variable(
                store_name + '_value',
                var_type=params['var_type'],
                dimensions=shape,
                chunksizes=chunksizes,
                simtk_unit=params['simtk_unit'],
            )

        value_store.initialize()

        # the value
        self.attribute_list[attribute] = value_store
        attribute_idx = self.storage.attributes.index[attribute.__uuid__]
        self.storage.attributes.vars['cache'][attribute_idx] = value_store

        # use the cache and function of the CV to fill the store when it is made
        if not allow_incomplete:

            indices = self.vars['uuid'][:]

            for pos, idx in enumerate(indices):

                proxy = LoaderProxy.new(key_store, idx)

                # value = attribute._cache_dict._get(proxy)
                #
                # if value is None:
                #     # not in cache so compute it if possible
                #     if attribute._eval_dict:
                #         value = attribute._eval_dict([proxy])[0]
                #     else:
                #         value = None

                value = attribute(proxy)

                if value is not None:
                    value_store.vars['value'][pos] = value
                    value_store.cache[pos] = value

        attribute.set_cache_store(value_store)
        return value_store
