"""
author
"""

import abc
import logging
import os.path
from collections import OrderedDict
from uuid import UUID

import netCDF4
import numpy as np
from .dictify import UUIDObjectJSON
from .stores import NamedObjectStore, ObjectStore, PseudoAttributeStore
from .proxy import LoaderProxy

import sys
if sys.version_info > (3, ):
    unicode = str

logger = logging.getLogger(__name__)
init_log = logging.getLogger('openpathsampling.initialization')

try:
    from simtk import unit as u
except ImportError:
    HAS_SIMTK_UNIT = False
else:
    HAS_SIMTK_UNIT = True


# ==============================================================================
# Extended NetCDF Storage for multiple forked trajectories
# ==============================================================================

class NetCDFPlus(netCDF4.Dataset):
    """
    Extension of the python netCDF wrapper for easier storage of python objects
    """
    support_simtk_unit = HAS_SIMTK_UNIT

    @property
    def _netcdfplus_version_(self):
        import openpathsampling.netcdfplus.version as v
        version = v.short_version
        return version

    _type_conversion = {
        'float': np.float32,
        'int': np.int32,
        'long': np.int64,
        'index': np.int32,
        'length': np.int32,
        'bool': np.int16,
        'str': str,
        'json': str,
        'jsonobj': str,
        'numpy.float32': np.float32,
        'numpy.float64': np.float64,
        'numpy.int8': np.int8,
        'numpy.int16': np.int16,
        'numpy.int32': np.int32,
        'numpy.int64': np.int64,
        'numpy.uint8': np.uint8,
        'numpy.uint16': np.uint16,
        'numpy.uint32': np.uint32,
        'numpy.uint64': np.uint64,
        'store': str,
        'obj': np.int32,
        'lazyobj': np.int32,
        'uuid': str
    }

    class ValueDelegate(object):
        """
        Value delegate for objects that implement __getitem__ and __setitem__

        It will basically just wrap values that are used in a dict like
        structure with getter and setter function to allow easier conversion

        delegate[x] is equivalent to delegate.getter(delegate.variable[x])

        Attributes
        ----------
        variable : dict-like
            the dict to be wrapped
        getter : function
            the function applied to results from running the __getitem__
            on the variable
        setter : function
            the function applied to the value to be stored using __setitem__
            on the variable
        store : openpathsampling.netcdfplus.ObjectStore
            a reference to an object store used for convenience in some cases

        """

        def __init__(self, variable, getter=None, setter=None, store=None):
            self.variable = variable
            self.store = store

            if setter is None:
                # None should not be used
                setter = lambda v: v
                self.__setitem__ = variable.__setitem__

            if getter is None:
                getter = lambda v: v
                self.__getitem__ = variable.__getitem__

            self.getter = getter
            self.setter = setter

            try:
                from simtk import unit as u
            except ImportError:
                self.support_simtk_unit = False

        def __setitem__(self, key, value):
            self.variable[key] = self.setter(value)

        def __getitem__(self, key):
            # print(self.variable[key])
            # print(type(self.variable[key]))
            return self.getter(self.variable[key])

        def __getattr__(self, item):
            return getattr(self.variable, item)

        def __str__(self):
            return str(self.variable)

        def __repr__(self):
            return repr(self.variable)

        def __len__(self):
            return len(self.variable)

    @property
    def objects(self):
        """
        Return a dictionary of all objects stored.

        This is similar to the netcdf `.variables` for all stored variables.
        This allows to write `storage.objects['samples'][idx]` like we
        write `storage.variables['ensemble_json'][idx]`
        """
        return self._stores

    def find_store(self, obj):
        """
        Return the default store used for an storable object

        Parameters
        ----------
        obj : :class:`openpathsampling.netcdfplus.StorableObject`
            the storable object to be tested

        Returns
        -------
        :class:`openpathsampling.netcdfplus.ObjectStore`
            the store that is used by default to store the given storable obj
        """

        if type(obj) is type or type(obj) is abc.ABCMeta:
            if obj not in self._obj_store:
                raise ValueError(
                    'Objects of class "%s" are not storable in this store.' %
                    obj.__name__)

            return self._obj_store[obj]
        else:
            if obj.__class__ not in self._obj_store:
                raise ValueError(
                    'Objects of class "%s" are not storable in this store.' %
                    obj.__class__.__name__)

            return self._obj_store[obj.__class__]

    def update_storable_classes(self):
        self.simplifier.update_class_list()

    def _create_storages(self):
        """
        Function to be called automatically to register all object stores

        This will usually only be called in subclassed storages.
        """
        # todo: add CVStore, rename to attribute
        pass

    def __init__(self, filename, mode=None, fallback=None):
        """
        Create a storage for complex objects in a netCDF file

        Parameters
        ----------
        filename : string
            filename of the netcdf file to be used or created
        mode : str
            the mode of file creation, one of 'w' (write), 'a' (append) or
            'r' (read-only) None, which will append any existing files
            (equal to append), is the default.
        fallback : :class:`openpathsampling.Storage`
            the _fallback_ storage to be loaded from if an object is not present
            in this storage. By default you will not try to resave objects
            that could be found in the fallback. Note that the fall back does
            only work if `use_uuid` is enabled

        Notes
        -----
        A single file can be opened by multiple storages, but only one can be
        used for writing

        """

        if mode is None:
            mode = 'a'

        self.mode = mode

        exists = os.path.isfile(filename)
        if exists and mode == 'a':
            logger.info(
                "Open existing netCDF file '%s' for appending - "
                "appending existing file", filename)
        elif exists and mode == 'w':
            logger.info(
                "Create new netCDF file '%s' for writing - "
                "deleting existing file", filename)
        elif not exists and mode == 'a':
            logger.info(
                "Create new netCDF file '%s' for appending - "
                "appending non-existing file", filename)
        elif not exists and mode == 'w':
            logger.info(
                "Create new netCDF file '%s' for writing - "
                "creating new file", filename)
        elif not exists and mode == 'r':
            logger.info(
                "Open existing netCDF file '%s' for reading - "
                "file does not exist", filename)
            raise RuntimeError("File '%s' does not exist." % filename)
        elif exists and mode == 'r':
            logger.info(
                "Open existing netCDF file '%s' for reading - "
                "reading from existing file", filename)

        self._filename = os.path.abspath(filename)
        self.fallback = fallback

        # this can be set to false to re-store objects present in the fallback
        self.exclude_from_fallback = True

        # this can be set to false to re-store proxies from other stores
        self.exclude_proxy_from_other = False

        # call netCDF4-python to create or open .nc file
        super(NetCDFPlus, self).__init__(filename, mode)

        self._setup_class()

        if mode == 'w':
            logger.info("Setup netCDF file and create variables")

            self.setncattr('format', 'netcdf+')
            self.setncattr('ncplus_version', self._netcdfplus_version_)

            self.write_meta()

            # add shared scalar dimension for everyone
            self.create_dimension('scalar', 1)
            self.create_dimension('pair', 2)

            self.setncattr('use_uuid', 'True')

            self._create_simplifier()

            # create the store that holds stores
            store_stores = NamedObjectStore(ObjectStore)
            store_stores.name = 'stores'
            self.register_store('stores', store_stores)
            self.stores.initialize()
            self.stores.set_caching(True)
            self.update_delegates()

            # now create all storages in subclasses
            self._create_storages()

            self.create_store('attributes', PseudoAttributeStore())

            # call the subclass specific initialization
            self._initialize()

            # this will create all variables in the storage for all new
            # added stores this is often already call inside of _initialize.
            # If not we just make sure
            self.finalize_stores()

            logger.info("Finished setting up netCDF file")

            self.sync()

        elif mode == 'a' or mode == 'r+' or mode == 'r':
            logger.debug("Restore the dict of units from the storage")

            self.check_version()

            # self.reference_by_uuid = hasattr(self, 'use_uuid')
            # self.reference_by_uuid = True

            self._create_simplifier()

            # open the store that contains all stores
            self.register_store('stores', NamedObjectStore(ObjectStore))
            self.stores.set_caching(True)
            self.create_variable_delegate('stores_json')
            self.create_variable_delegate('stores_name')
            self.create_variable_delegate('stores_uuid')

            self.stores.restore()

            # Create a dict of simtk.Unit() instances for all netCDF.Variable()
            for variable_name in self.variables:
                variable = self.variables[variable_name]

                if self.support_simtk_unit:
                    if hasattr(variable, 'unit_simtk'):
                        unit_dict = self.simplifier.from_json(
                            getattr(variable, 'unit_simtk'))
                        if unit_dict is not None:
                            unit = self.simplifier.unit_from_dict(unit_dict)
                        else:
                            unit = self.simplifier.unit_from_dict(u.Unit({}))

                        self.units[str(variable_name)] = unit

            # register all stores that are listed in self.stores

            for store in self.stores:
                if store is not None:
                    logger.debug("Register store %s in the storage" % store.name)
                    self.register_store(store.name, store)
                    store.register(self, store.name)

            self.update_delegates()
            self._restore_storages()

            # only if we have a new style file
            if hasattr(self, 'attributes'):
                for attribute, store in zip(
                        self.attributes,
                        self.attributes.vars['cache']
                ):
                    if store is not None:
                        key_store = self.attributes.key_store(attribute)
                        key_store.attribute_list[attribute] = store

            # call the subclass specific restore in case there is more stuff
            # to prepare
            self._restore()

        self.set_auto_mask(False)
        # self.set_always_mask(False)  ## didn't fix; errors older versions



    def _create_simplifier(self):
        self.simplifier = UUIDObjectJSON(self)

    @property
    def filename(self):
        return self._filename

    @property
    def file_size(self):
        return os.path.getsize(self.filename)

    @property
    def file_size_str(self):
        current = float(self.file_size)
        output_prefix = ''
        for prefix in ["k", "M", "G"]:
            if current >= 1024:
                output_prefix = prefix
                current /= 1024.0
        return "{0:.2f}{1}B".format(current, output_prefix)

    @staticmethod
    def _cmp_version(v1, v2):
        # we only look at x.y.z parts
        def version_parts(v):
            return v.split('-')[0].split('+')[0].split('.')[:3]
        q1 = version_parts(v1)
        q2 = version_parts(v2)
        # q1 = v1.split('.')[:3]
        # q2 = v2.split('.')[:3]
        # q1 = v1.split('-')[0].split('.')
        # q2 = v2.split('-')[0].split('.')
        for v1, v2 in zip(q1, q2):
            if int(v1) > int(v2):
                return +1
            elif int(v1) < int(v2):
                return -1

        return 0

    def check_version(self):
        try:
            s1 = self.getncattr('ncplus_version')
        except AttributeError:
            logger.info(
                'Using netcdfplus Pre 1.0 version. '
                'No version detected using 0.0.0')
            s1 = '0.0.0'

        s2 = self._netcdfplus_version_

        cp = self._cmp_version(s1, s2)

        if cp != 0:
            logger.info('Loading different netcdf version. '
                        'Installed version is '
                        '%s and loaded version is %s' % (s2, s1))
            if cp > 0:
                logger.info(
                    'Loaded version is newer consider upgrading your '
                    'conda package!')
            else:
                logger.info(
                    'Loaded version is older. Should be no problem other then '
                    'missing features and information')

    def write_meta(self):
        pass

    def _setup_class(self):
        """
        Sets the basic properties for the storage
        """
        self._stores = OrderedDict()
        self._objects = {}
        self._obj_store = {}
        self._storages_base_cls = {}
        self.vars = dict()
        self.units = dict()

    def create_store(self, name, store, register_attr=True):
        """
        Create a special variable type `obj.name` that can hold storable objects

        Parameters
        ----------
        name : str
            the name of the store inside this storage
        store : :class:`openpathsampling.netcdf.ObjectStore`
            the store to be added to this storage
        register_attr : bool
            if `True` the store will be added to the storage as an
             attribute with name `name`

        """
        self.register_store(name, store, register_attr=register_attr)
        store.name = name
        self.stores.save(store)

    def finalize_stores(self):
        """
        Run initializations for all added stores.

        This will make sure that all previously added stores are now useable.
        If you add more stores you need to call this again. The reason this is
        done at all is that stores might reference each other and so no
        unique order of creation can be found. Thus you first create stores
        with all their dependencies and then finalize all of them together.
        """
        for store in self._stores.values():
            if not store.is_created():
                logger.info("Initializing store '%s'" % store.name)
                store.initialize()

        for store in self._stores.values():
            if not store.is_created():
                logger.info("Initializing store '%s'" % store.name)
                store.initialize()

        self.update_delegates()
        self.simplifier.update_class_list()

    def register_store(self, name, store, register_attr=True):
        """
        Add a object store to the file

        An object store is a special type of variable that allows to store
        python objects

        Parameters
        ----------
        name : str
            the name of the store under which the objects are accessible
            like `store.{name}`
        store : :class:`openpathsampling.storages.ObjectStore`
            instance of the object store
        register_attr : bool
            if set to false the store will not be accesible as an attribute.
            `True` is the default.
        """
        store.register(self, name)

        if register_attr:
            if hasattr(self, name):
                raise ValueError('Store name %s is already in use!' % name)

            setattr(self, store.prefix, store)

        self._stores[name] = store

        if store.content_class is not None:
            self._objects[store.content_class] = store

            self._obj_store[store.content_class] = store
            self._obj_store.update(
                {cls: store for cls in store.content_class.descendants()})

    def _initialize(self):
        """
        Function run after a new file is created.

        This is used to setup all variables in the storage
        """
        pass

    def _restore(self):
        """
        Function run after an existing file is opened.

        This is used in special storages to complete reading existing files.
        """
        pass

    def __repr__(self):
        return "Storage @ '" + self.filename + "'"

    def __getattr__(self, item):
        try:
            return self.__dict__[item]
        except KeyError:
            try:
                return self.__class__.__dict__[item]
            except KeyError:
                raise AttributeError(
                    "'{cls}' object has no attribute '{itm}'".format(
                        cls=self.__class__, itm=item
                    )
                )

    def __setattr__(self, key, value):
        self.__dict__[key] = value

    def _init_storages(self):
        """
        Run the initialization on all added classes

        Notes
        -----
        Only runs when the storage is created.
        """

        for storage in self._stores.values():
            storage.initialize()

        self.update_delegates()

    def _restore_storages(self):
        """
        Run the restore method on all added classes

        Notes
        -----
        Only runs when an existing storage is opened.
        """

        for storage in self._stores.values():
            storage.restore()
            storage._created = True

    def list_stores(self):
        """
        Return a list of registered stores

        Returns
        -------
        list of str
            list of stores that can be accessed using `storage.[store]`
        """
        return [store.prefix for store in self._stores.values()]

    def list_storable_objects(self):
        """
        Return a list of storable object base classes

        Returns
        -------
        list of type
            list of base classes that can be stored using `storage.save(obj)`
        """
        return [
            store.content_class
            for store in self.objects.values()
            if store.content_class is not None]

    def save(self, obj, idx=None):
        """
        Save a storable object into the correct Storage in the netCDF file

        Parameters
        ----------
        obj : :class:`StorableObject`
            the object to store
        idx : str
            the name to be stored by

        Returns
        -------
        str
            the class name of the BaseClass of the stored object, which is
            needed when loading the object to identify the correct storage
        """

        if type(obj) is list:
            # a list of objects will be stored one by one
            return [self.save(part, idx) for part in obj]

        elif type(obj) is tuple:
            # a tuple will store all parts
            return [self.save(part, idx) for part in obj]

        elif obj.__class__ in self._obj_store:
            # to store we just check if the base_class is present in the
            # storages also we assume that if a class has no base_cls
            store = self.find_store(obj)
            store_idx = self.stores.index[store.__uuid__]
            return store, store_idx, store.save(obj, idx)

        # Could not save this object.
        raise RuntimeWarning("Objects of type '%s' cannot be stored!" %
                             obj.__class__.__name__)

    def __contains__(self, item):
        if type(item) is list:
            # a list of objects will be stored one by one
            return [part in self for part in item]

        elif type(item) is tuple:
            # a tuple will store all parts
            return tuple([part in self for part in item])

        elif item.__class__ in self._obj_store:
            # to store we just check if the base_class is present in the
            # storages also we assume that if a class has no base_cls
            store = self.find_store(item)
            return item in store

        return False

    def load(self, uuid):
        """
        Load an object from the storage

        Parameters
        ----------
        uuid : uuid.UUID
            the uuid to be loaded

        Returns
        -------
        :class:`openpathsampling.netcdfplus.StorableObject`
            the object loaded from the storage

        Notes
        -----
        this only works in storages with uuids otherwise load directly from the
        substores
        """

        for store in self.objects.values():
            if uuid in store.index:
                return store[uuid]

        raise KeyError("UUID %s not found in storage" % uuid)

    def idx(self, obj):
        """
        Return the index used to store the given object in this storage

        Parameters
        ----------
        obj : object
            The stored object from which the index is to be returned
        """
        if hasattr(obj, 'base_cls'):
            store = self._objects[obj.base_cls]
            return store.idx(obj)

    def repr_json(self, obj):
        """
        Return the JSON representation in the storage if available

        Parameters
        ----------
        obj : :class:`openpathsampling.netcdfplus.StorableObject`

        Returns
        -------
        str
            the JSON string (usually in unicode) from the storage
        """
        if hasattr(obj, 'base_cls'):
            store = self._objects[obj.base_cls]

            if store.json:
                return store.variables['json'][store.idx(obj)]

        return None

    def create_dimension(self, dim_name, size=None):
        """
        Initialize a new dimension in the storage.
        Wraps the netCDF createDimension

        Parameters
        ----------
        dim_name : str
            the name for the new dimension
        size : int
            the number of elements in this dimension. None (default) means
            an infinite dimension that extends when more objects are stored

        """
        if dim_name not in self.dimensions:
            self.createDimension(dim_name, size)

    def cache_image(self):
        """
        Return an dict containing information about all caches

        Returns
        -------
        dict
            a nested dict containing information about the number and types of
            cached objects
        """
        image = {
            'weak': {},
            'strong': {},
            'total': {},
            'file': {},
            'index': {}
        }

        total_strong = 0
        total_weak = 0
        total_file = 0
        total_index = 0

        for name, store in self.objects.items():
            size = store.cache.size
            count = store.cache.count
            profile = {
                'count': count[0] + count[1],
                'count_strong': count[0],
                'count_weak': count[1],
                'max': size[0],
                'size_strong': size[0],
                'size_weak': size[1],
            }
            total_strong += count[0]
            total_weak += count[1]
            total_file += len(store)
            total_index += len(store.index)
            image[name] = profile
            image['strong'][name] = count[0]
            image['weak'][name] = count[1]
            image['total'][name] = count[0] + count[1]
            image['file'][name] = len(store)
            #            if hasattr(store, 'index'):
            image['index'][name] = len(store.index)
        # else:
        #                image['index'][name] = 0

        image['full'] = total_weak + total_strong
        image['total_strong'] = total_strong
        image['total_weak'] = total_weak
        image['file'] = total_file
        image['index'] = total_index

        return image

    def get_var_types(self):
        """
        List all allowed variable type to be used in `create_variable`

        Returns
        -------
        list of str
            the list of variable types
        """
        types = list(NetCDFPlus._type_conversion.keys())
        types += ['obj.' + x for x in self.objects.keys()]
        types += ['lazyobj.' + x for x in self.objects.keys()]
        types += ['uuid.' + x for x in self.objects.keys()]
        types += ['lazyuuid.' + x for x in self.objects.keys()]
        return sorted(types)

    @staticmethod
    def var_type_to_nc_type(var_type):
        """
        Return the compatible netCDF variable type for var_type

        Parameters
        ----------
        var_type :

        Returns
        -------
        object
            A object of netcdf compatible varible types
        """
        if 'obj.' in var_type:
            nc_type = str
        elif 'uuid.' in var_type:
            nc_type = str
        else:
            nc_type = NetCDFPlus._type_conversion[var_type]

        return nc_type

    def create_type_delegate(self, var_type):
        """
        Create a variable value delegator for var_type

        The delegator will convert automatically between the given variable type
        and the netcdf compatible one

        Parameters
        ----------
        var_type : str
            the variable type

        Returns
        -------
        NetCDFPlus.Value_Delegate
            the delegator instance
        """
        getter = None
        setter = None
        store = None

        if 'obj.' in var_type or 'uuid.' in var_type:
            store_name = str(var_type.split('.')[1])
            store = self._stores[store_name]

            base_type = store.content_class

            # get_is_iterable = lambda v: \
            #     v.base_cls is not base_type if hasattr(v, 'base_cls') else \
            #     hasattr(v, '__iter__')

            get_numpy_iterable = lambda v: isinstance(v, np.ndarray)

            set_is_iterable = lambda v: \
                v.base_cls is not base_type if hasattr(v, 'base_cls') else \
                hasattr(v, '__iter__')

        if var_type == 'int':
            getter = lambda v: v.tolist()
            setter = lambda v: np.array(v)

        elif var_type == 'bool':
            getter = lambda v: v.astype(bool).tolist()
            setter = lambda v: np.array(v, dtype=np.int8)

        elif var_type == 'index':
            getter = lambda v: \
                [None if int(w) < 0 else int(w) for w in v.tolist()] \
                if hasattr(v, '__iter__') else None if int(v) < 0 else int(v)
            setter = lambda v: \
                [-1 if w is None else w for w in v] \
                if hasattr(v, '__iter__') else -1 if v is None else v

        elif var_type == 'float':
            getter = lambda v: v.tolist()
            setter = lambda v: np.array(v)

        elif var_type.startswith('numpy.'):
            pass

        elif var_type == 'jsonobj':
            setter = lambda v: self.simplifier.to_json_object(v)
            getter = lambda v: self.simplifier.from_json(v)

        elif var_type == 'json':
            setter = lambda v: self.simplifier.to_json(v)
            getter = lambda v: self.simplifier.from_json(v)

        elif var_type.startswith('obj.'):
            getter = lambda v: [
                None if w[0] == '-' else store.load(int(UUID(w)))
                for w in v
            ] if get_numpy_iterable(v) else \
                None if v[0] == '-' else store.load(int(UUID(v)))

            setter = lambda v: \
                ''.join(['-' * 36 if w is None else str(UUID(int=store.save(w)))
                         for w in list.__iter__(v)]) \
                if set_is_iterable(v) else \
                '-' * 36 if v is None else str(UUID(int=store.save(v)))

        elif var_type.startswith('lazyobj.'):
            getter = lambda v: [
                None if w[0] == '-' else LoaderProxy.new(store, int(UUID(w)))
                for w in v
            ] if isinstance(v, np.ndarray) else \
                None if v[0] == '-' else LoaderProxy.new(store, int(UUID(v)))

            setter = lambda v: \
                ''.join([
                    '-' * 36 if w is None else str(UUID(int=store.save(w)))
                    for w in list.__iter__(v)
                ]) if set_is_iterable(v) else \
                '-' * 36 if v is None else str(UUID(int=store.save(v)))

        elif var_type == 'uuid':
            getter = lambda v: \
                [None if w[0] == '-' else int(UUID(w)) for w in v] \
                if type(v) is not unicode else None \
                if v[0] == '-' else int(UUID(v))
            setter = lambda v: \
                ''.join([
                    '-' * 36 if w is None else str(UUID(int=w))
                    for w in list.__iter__(v)
                ]) if hasattr(v, '__iter__') else \
                '-' * 36 if v is None else str(UUID(int=v))

        elif var_type == 'store':
            setter = lambda v: v.prefix
            getter = lambda v: self.objects[v]

        return getter, setter, store

    to_uuid_chunks = staticmethod(
        lambda x: [x[i:i + 36] for i in range(0, len(x), 36)])

    def create_variable_delegate(self, var_name):
        """
        Create a delegate property that wraps the netcdf.Variable and takes care
        of type conversions

        Parameters
        ----------
        var_name : str
            the name of the variable for which a delegate should be created

        """

        if var_name not in self.vars:
            var = self.variables[var_name]

            if not hasattr(var, 'var_type'):
                return

            getter, setter, store = self.create_type_delegate(var.var_type)

            to_uuid_chunks = NetCDFPlus.to_uuid_chunks
            # to_uuid_chunks34 = NetCDFPlus.to_uuid_chunks34

            if hasattr(var, 'var_vlen'):
                if var.var_type.startswith('obj.'):
                    getter = lambda v: [[
                        None if u[0] == '-' else store.load(int(UUID(u)))
                        for u in to_uuid_chunks(w)
                        ] for w in v
                    ] if isinstance(v, np.ndarray) else [
                        None if u[0] == '-' else store.load(int(UUID(u)))
                        for u in to_uuid_chunks(v)
                    ]
                elif var.var_type.startswith('lazyobj.'):
                    getter = lambda v: [[
                        None if u[0] == '-' else
                        LoaderProxy.new(store, int(UUID(u)))
                        for u in to_uuid_chunks(w)] for w in v
                    ] if isinstance(v, np.ndarray) else [
                        None if u[0] == '-' else
                        LoaderProxy.new(store, int(UUID(u)))
                        for u in to_uuid_chunks(v)
                    ]

            if True or self.support_simtk_unit:
                if hasattr(var, 'unit_simtk'):
                    if var_name not in self.units:
                        self.update_simtk_unit(var_name)

                    unit = self.units[var_name]

                    def _get(my_getter):
                        import simtk.unit as u
                        if my_getter is None:
                            return lambda v: u.Quantity(v, unit)
                        else:
                            return lambda v: u.Quantity(my_getter(v), unit)

                    def _set(my_setter):
                        if my_setter is None:
                            return lambda v: v / unit
                        else:
                            return lambda v: my_setter(v / unit)

                    getter = _get(getter)
                    setter = _set(setter)

            if True:
                if hasattr(var, 'maskable'):
                    def _get2(my_getter):
                        return lambda v: [
                            None if hasattr(w, 'mask') else
                            my_getter(w) for w in v
                        ] if type(v) is not str and len(v.shape) > 0 else \
                            None if hasattr(v, 'mask') else my_getter(v)

                    if getter is not None:
                        getter = _get2(getter)
                    else:
                        getter = _get2(lambda v: v)

            delegate = NetCDFPlus.ValueDelegate(var, getter, setter, store)

            # this is a trick to speed up the s/getter. If we do not need
            # to _cast_ because of python objects of units we can copy
            # the s/getter of the original var which is still bound to the
            # right object

            self.vars[var_name] = delegate

        else:
            raise ValueError("Variable '%s' is already taken!" % var_name)

    def create_variable(self, var_name,
                        var_type,
                        dimensions,
                        description=None,
                        chunksizes=None,
                        simtk_unit=None,
                        maskable=False):
        """
        Create a new variable in the netCDF storage.

        This is just a helper function to structure the code better and add
        some convenience to creating more complex variables

        Parameters
        ==========
        var_name : str
            The name of the variable to be created
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

        ncfile = self

        if type(dimensions) is str:
            dimensions = [dimensions]

        dimensions = list(dimensions)

        new_dimensions = dict()
        for ix, dim in enumerate(dimensions):
            if type(dim) is int:
                dimensions[ix] = var_name + '_dim_' + str(ix)
                new_dimensions[dimensions[ix]] = dim

        if dimensions[-1] == '...':
            # last dimension is simply [] so we allow arbitrary length
            # and remove the last dimension
            variable_length = True
            dimensions = dimensions[:-1]
        else:
            variable_length = False

        if var_type == 'obj' or var_type == 'lazyobj':
            dimensions.append('pair')
            if chunksizes is not None:
                chunksizes = tuple(list(chunksizes) + [2])

        nc_type = self.var_type_to_nc_type(var_type)

        for dim_name, size in new_dimensions.items():
            ncfile.create_dimension(dim_name, size)

        dimensions = tuple(dimensions)

        # if chunk sizes are strings then replace it by
        # the actual size of the dimension
        if chunksizes is not None:
            chunksizes = list(chunksizes)
            for ix, dim in enumerate(chunksizes):
                if dim == -1:
                    chunksizes[ix] = len(ncfile.dimensions[dimensions[ix]])

                if type(dim) is str:
                    chunksizes[ix] = len(ncfile.dimensions[dim])

            chunksizes = tuple(chunksizes)

        if variable_length:
            vlen_t = ncfile.createVLType(nc_type, var_name + '_vlen')
            ncvar = ncfile.createVariable(
                var_name, vlen_t, dimensions, chunksizes=chunksizes
            )

            setattr(ncvar, 'var_vlen', 'True')
        else:
            ncvar = ncfile.createVariable(
                var_name, nc_type, dimensions, chunksizes=chunksizes,
            )

        setattr(ncvar, 'var_type', var_type)

        if self.support_simtk_unit and simtk_unit is not None:

            import simtk.unit as u

            if isinstance(simtk_unit, u.Unit):
                unit_instance = simtk_unit
                symbol = unit_instance.get_symbol()
            elif isinstance(simtk_unit, u.BaseUnit):
                unit_instance = u.Unit({simtk_unit: 1.0})
                symbol = unit_instance.get_symbol()
            elif type(simtk_unit) is str and hasattr(u, simtk_unit):
                unit_instance = getattr(u, simtk_unit)
                symbol = unit_instance.get_symbol()
            else:
                raise NotImplementedError(
                    'Unit by abbreviated string representation '
                    'is not yet supported')

            json_unit = self.simplifier.unit_to_json(unit_instance)

            # store the unit in the dict inside the Storage object
            self.units[var_name] = unit_instance

            # Define units for a float variable
            setattr(ncvar, 'unit_simtk', json_unit)
            setattr(ncvar, 'unit', symbol)

        if maskable:
            setattr(ncvar, 'maskable', 'True')

        if description is not None:
            if type(dimensions) is str:
                dim_names = [dimensions]
            else:
                dim_names = ['#ix{0}:{1}'.format(*p) for p in
                             enumerate(dimensions)]

            idx_desc = '[' + ']['.join(dim_names) + ']'
            description = var_name + idx_desc + ' is ' + \
                description.format(idx=dim_names[0], ix=dim_names)

            # Define long (human-readable) names for variables.
            setattr(ncvar, "long_str", description)

        self.update_delegates()

        return ncvar

    def update_delegates(self):
        """
        Updates the set of delegates in `self.vars`

        Should be called after new variables have been created or loaded.
        """
        for name in self.variables:
            if name not in self.vars:
                self.create_variable_delegate(name)

    @staticmethod
    def get_value_parameters(value):
        """
        Compute netcdfplus compatible parameters to store a value

        Parameters
        ----------
        value

        Returns
        -------
        dict
            A dictionary containing the approriate input parameters for
            `var_type`, `dimensions`, `simtk_unit`

        Notes
        -----
        This is a utility function to create a CV using a template

        """

        dimensions = None
        storable = True
        simtk_unit = None

        test_value = value
        test_type = value

        if NetCDFPlus.support_simtk_unit:
            import simtk.unit as u

            if type(test_type) is u.Quantity:
                # could be a Quantity([..])
                simtk_unit = test_type.unit
                test_type = test_type._value
        else:
            u = None

        if type(test_type) is np.ndarray:
            dimensions = test_type.shape
        else:
            if hasattr(test_value, '__len__'):
                dimensions = len(test_value)
                test_type = test_value[0]
                if NetCDFPlus.support_simtk_unit and type(test_type) \
                        is u.Quantity:
                    for val in test_value:
                        if isinstance(val._value, type(test_value._value)):
                            # all values must be of same type
                            storable = False
                else:
                    for val in test_value:
                        if type(val) is not type(test_value):
                            # all values must be of same type
                            storable = False

            if NetCDFPlus.support_simtk_unit and type(test_type) is u.Quantity:
                # could also be [Quantity, ...]
                simtk_unit = test_type.unit
                test_type = test_type._value

        if storable:
            var_type = NetCDFPlus.identify_var_type(test_type)
            return {
                'var_type': var_type,
                'dimensions': dimensions,
                'simtk_unit': simtk_unit
            }

        return {
        }

    @staticmethod
    def identify_var_type(instance):
        """
        Identify common python and numpy types

        Parameters
        ----------
        instance
            python variable instance to be tested for it numeric type

        Returns
        -------
        str
            a string representation of the variable type

        """
        ty = type(instance)

        known_types = [float, int, bool, str]

        if ty in known_types:
            return ty.__name__
        elif hasattr(instance, 'dtype'):
            return 'numpy.' + instance.dtype.type.__name__
        else:
            return 'None'
