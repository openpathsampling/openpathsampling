"""
author
"""

import logging

from dictify import StorableObjectJSON
from proxy import LoaderProxy

from objects import NamedObjectStore, ObjectStore

import numpy as np
import netCDF4
import os.path
import abc

logger = logging.getLogger(__name__)
init_log = logging.getLogger('openpathsampling.initialization')


# =============================================================================================
# Extended NetCDF Storage for multiple forked trajectories
# =============================================================================================

class NetCDFPlus(netCDF4.Dataset):
    """
    Extension of the python netCDF wrapper for easier storage of python objects
    """
    support_simtk_unit = True

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
        'lazyobj': np.int32
    }

    class ValueDelegate(object):
        """
        Value delegate for objects that implement __getitem__ and __setitem__

        It will basically just wrap values that are used in a dict like structure
        with getter and setter function to allow easier conversion

        delegate[x] is equivalent to delegate.getter(delegate.variable[x])

        Attributes
        ----------
        variable : dict-like
            the dict to be wrapped
        getter : function
            the function applied to results from running the __getitem__ on the variable
        setter : function
            the function applied to the value to be stored using __setitem__ on the variable
        store : openpathsampling.netcdfplus.ObjectStore
            a reference to an object store used for convenience in some cases

        """

        def __init__(self, variable, getter=None, setter=None, store=None):
            self.variable = variable
            self.store = store

            if setter is None:
                setter = lambda v: v
            self.setter = setter

            if getter is None:
                getter = lambda v: v
            self.getter = getter

        def __setitem__(self, key, value):
            self.variable[key] = self.setter(value)

        def __getitem__(self, key):
            return self.getter(self.variable[key])

        def __getattr__(self, item):
            return getattr(self.variable, item)

        def __str__(self):
            return str(self.variable)

        def __repr__(self):
            return repr(self.variable)

    class KeyDelegate(object):
        """
        Value delegate for objects that implement __getitem__ and __setitem__

        It will basically just wrap keys for objects that are used in a dict like structure
        with getter and setter function to allow easier conversion

        delegate[x] is equivalent to delegate[x.idx(store)]

        Attributes
        ----------
        variable : dict-like
            the dict to be wrapped
        store : openpathsampling.netcdfplus.ObjectStore
            a reference to an object store used

        """

        def __init__(self, variable, store):
            self.variable = variable
            self.store = store

        def __setitem__(self, key, value):
            if hasattr(key, '__iter__'):
                idxs = [item if type(item) is int else self.store.index[item] for item in key]
                sorted_idxs = list(set(idxs))
                sorted_values = [value[idxs.index(val)] for val in sorted_idxs]
                self.variable[sorted_idxs] = sorted_values

            else:
                self.variable[key if type(key) is int else self.store.index[key]] = value

        def __getitem__(self, key):
            if hasattr(key, '__iter__'):
                idxs = [item if type(item) is int else self.store.index[item] for item in key]
                sorted_idxs = sorted(list(set(idxs)))

                sorted_values = self.variable[sorted_idxs]
                return [sorted_values[sorted_idxs.index(idx)] for idx in idxs]
            else:
                return self.variable[key if type(key) is int else self.store.index[key]]

    @property
    def objects(self):
        """
        Return a dictionary of all objects stored.

        This is similar to the netcdf `.variables` for all stored variables. This
        allows to write `storage.objects['samples'][idx]` like we
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
                raise ValueError('Objects of class "%s" are not storable in this store.' % obj.__name__)

            return self._obj_store[obj]
        else:
            if obj.__class__ not in self._obj_store:
                raise ValueError('Objects of class "%s" are not storable in this store.' % obj.__class__.__name__)

            return self._obj_store[obj.__class__]

    def update_storable_classes(self):
        self.simplifier.update_class_list()

    def _create_storages(self):
        """
        Function to be called automatically to register all object stores

        This will usually only be called in subclassed storages.
        """
        pass

    def __init__(self, filename, mode=None):
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

        Notes
        -----
        A single file can be opened by multiple storages, but only one can be used for writing

        """

        if mode is None:
            mode = 'a'

        exists = os.path.isfile(filename)
        if exists and mode == 'a':
            logger.info("Open existing netCDF file '%s' for appending - appending existing file", filename)
        elif exists and mode == 'w':
            logger.info("Create new netCDF file '%s' for writing - deleting existing file", filename)
        elif not exists and mode == 'a':
            logger.info("Create new netCDF file '%s' for appending - appending non-existing file", filename)
        elif not exists and mode == 'w':
            logger.info("Create new netCDF file '%s' for writing - creating new file", filename)
        elif not exists and mode == 'r':
            logger.info("Open existing netCDF file '%s' for reading - file does not exist", filename)
            raise RuntimeError("File '%s' does not exist." % filename)
        elif exists and mode == 'r':
            logger.info("Open existing netCDF file '%s' for reading - reading from existing file", filename)

        self.filename = filename

        # call netCDF4-python to create or open .nc file
        super(NetCDFPlus, self).__init__(filename, mode)

        self._setup_class()

        if mode == 'w':
            logger.info("Setup netCDF file and create variables")

            # add shared scalar dimension for everyone
            self.create_dimension('scalar', 1)
            self.create_dimension('pair', 2)

            # create the store that holds stores
            self.register_store('stores', NamedObjectStore(ObjectStore))
            self.stores._init()
            self.stores.set_caching(True)
            self.update_delegates()

            # now create all storages in subclasses
            self._create_storages()

            # call the subclass specific initialization
            self._initialize()

            # this will create all variables in the storage for all new added stores
            # this is often already call inside of _initialize. If not we just make sure
            self.finalize_stores()

            logger.info("Finished setting up netCDF file")

        elif mode == 'a' or mode == 'r+' or mode == 'r':
            logger.debug("Restore the dict of units from the storage")

            # open the store that contains all stores
            self.register_store('stores', NamedObjectStore(ObjectStore))
            self.stores.set_caching(True)
            self.create_variable_delegate('stores_json')
            self.create_variable_delegate('stores_name')

            # Create a dict of simtk.Unit() instances for all netCDF.Variable()
            for variable_name in self.variables:
                variable = self.variables[variable_name]

                if self.support_simtk_unit:
                    import simtk.unit as u
                    if hasattr(variable, 'unit_simtk'):
                        unit_dict = self.simplifier.from_json(getattr(variable, 'unit_simtk'))
                        if unit_dict is not None:
                            unit = self.simplifier.unit_from_dict(unit_dict)
                        else:
                            unit = self.simplifier.unit_from_dict(u.Unit({}))

                        self.units[str(variable_name)] = unit

            # register all stores that are listed in self.stores
            for store in self.stores:
                self.register_store(store.name, store)
                store.register(self, store.name)

            self.update_delegates()
            self._restore_storages()

            # call the subclass specific restore in case there is more stuff the prepare
            self._restore()

        self.sync()

    def _setup_class(self):
        """
        Sets the basic properties for the storage
        """
        self._stores = {}
        self._objects = {}
        self._obj_store = {}
        self._storages_base_cls = {}
        self.simplifier = StorableObjectJSON(self)
        self.vars = dict()
        self.units = dict()

    def create_store(self, name, store):
        """
        Create a special variable type `obj.name` that can hold storable objects

        Parameters
        ----------
        name : str
            the name of the store inside this storage
        store : :class:`openpathsampling.netcdf.ObjectStore`
            the store to be added to this storage

        """
        self.register_store(name, store)
        store.name = name
        self.stores.save(store)

    def finalize_stores(self):
        """
        Run initializations for all added stores.

        This will make sure that all previously added stores are now useable. If you add more stores you need to
        call this again. The reason this is done at all is that stores might reference each other and so no
        unique order of creation can be found. Thus you first create stores with all their dependencies and then
        finalize all of them together.
        """
        for store in self._stores.values():
            if not store._created:
                logger.info("Initializing store '%s'" % store.name)
                store._init()

        for store in self._stores.values():
            if not store._created:
                logger.info("Initializing store '%s'" % store.name)
                store._init()

        self.update_delegates()
        self.simplifier.update_class_list()

    def register_store(self, name, store, register_attr=True):
        """
        Add a object store to the file

        An object store is a special type of variable that allows to store python objects

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
                raise ValueError('Attribute name %s is already in use!' % name)

            setattr(self, store.prefix, store)

        self._stores[name] = store

        if store.content_class is not None:
            self._objects[store.content_class] = store

            self._obj_store[store.content_class] = store
            self._obj_store.update({cls: store for cls in store.content_class.descendants()})

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
        return self.__dict__[item]

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
            storage._init()

        self.update_delegates()

    def _restore_storages(self):
        """
        Run the restore method on all added classes

        Notes
        -----
        Only runs when an existing storage is opened.
        """

        for storage in self._stores.values():
            storage._restore()

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
        return [store.content_class for store in self.objects.values() if store.content_class is not None]

    def save(self, obj, idx=None):
        """
        Save a storable object into the correct Storage in the netCDF file

        Parameters
        ----------
        obj : the object to store

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
            # to store we just check if the base_class is present in the storages
            # also we assume that if a class has no base_cls
            store = self.find_store(obj)
            store_idx = self.stores.index[store]
            return store, store_idx, store.save(obj, idx)


        # Could not save this object.
        raise RuntimeWarning("Objects of type '%s' cannot be stored!" %
                             obj.__class__.__name__)

    def load(self, obj_type, idx):
        """
        Load an object of the specified type from the storage

        Parameters
        ----------
        obj_type : str or class
            the store or class of the base object to be loaded.

        Returns
        -------
        :class:`openpathsampling.netcdfplus.StorableObject`
            the object loaded from the storage

        Notes
        -----
        If you want to load a sub-classed Ensemble you need to load using
        `Ensemble` or `"Ensemble"` and not use the subclass
        """

        if obj_type in self._objects:
            store = self._objects[obj_type]
            return store.load(idx)
        elif obj_type in self._obj_store:
            # check if a store for the base_cls exists and use this one
            store = self._obj_store[obj_type]
            return store.load(idx)
        elif obj_type in self.simplifier.class_list:
            store = self._obj_store[self.simplifier.class_list[obj_type]]
            return store.load(idx)

        raise RuntimeError("No store registered to load variable type '%s'" % obj_type)

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
        if hasattr(obj, 'base_cls'):
            store = self._objects[obj.base_cls]

            if store.json:
                return store.variables['json'][store.idx(obj)]

        return None

    def clone_store(self, store_to_copy, new_storage):
        """
        Clone a store from one storage to another. Mainly used as a helper
        for the cloning of a store

        Parameters
        ----------
        store_to_copy : [..]Store
            the store to be copied
        new_storage : Storage
            the new Storage object

        """
        if type(store_to_copy) is str:
            storage_name = store_to_copy
        else:
            storage_name = store_to_copy.prefix

        copied_storages = 0

        for variable in self.variables.keys():
            if variable.startswith(storage_name + '_'):
                copied_storages += 1
                if variable not in new_storage.variables:
                    # collectivevariables have additional variables in the storage that need to be copied
                    # TODO: copy chunksizes?
                    var = self.variables[variable]
                    new_storage.createVariable(
                        variable,
                        str(var.dtype),
                        var.dimensions,
                        chunksizes=var.chunk
                    )
                    for attr in self.variables[variable].ncattrs():
                        setattr(
                            new_storage.variables[variable],
                            attr,
                            getattr(self.variables[variable], attr)
                        )

                    new_storage.variables[variable][:] = self.variables[variable][:]
                else:
                    for idx in range(0, len(self.variables[variable])):
                        new_storage.variables[variable][idx] = self.variables[variable][idx]

        if copied_storages == 0:
            raise RuntimeWarning(
                'Potential error in storage name. No storage variables copied from ' +
                storage_name
            )

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

        for name, store in self.objects.iteritems():
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
        types = NetCDFPlus._type_conversion.keys()
        types += ['obj.' + x for x in self.objects.keys()]
        types += ['lazyobj.' + x for x in self.objects.keys()]
        return sorted(types)

    @staticmethod
    def var_type_to_nc_type(var_type):
        """
        Return the compatible netCDF variable type for var_type

        Returns
        -------
        object
            A object of netcdf compatible varible types
        """
        if var_type.startswith('obj.') or var_type.startswith('lazyobj.'):
            nc_type = np.int32
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

        if var_type.startswith('obj.') or var_type.startswith('lazyobj.'):
            store = getattr(self, var_type.split('.')[1])
            base_type = store.content_class

            get_is_iterable = lambda v: \
                not v.base_cls is base_type if hasattr(v, 'base_cls') else hasattr(v, '__iter__')

            set_is_iterable = lambda v: \
                not v.base_cls is base_type if hasattr(v, 'base_cls') else hasattr(v, '__iter__')


        if var_type == 'int':
            getter = lambda v: v.tolist()
            setter = lambda v: np.array(v)

        elif var_type == 'bool':
            getter = lambda v: v.astype(np.bool).tolist()
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
            getter = lambda v: \
                [None if int(w) < 0 else store.load(int(w)) for w in v.tolist()] \
                    if get_is_iterable(v) else None if int(v) < 0 else store.load(int(v))
            setter = lambda v: \
                np.array([-1 if w is None else store.save(w) for w in v], dtype=np.int32) \
                    if set_is_iterable(v) else -1 if v is None else store.save(v)

        elif var_type == 'obj':
            # arbitrary object

            set_iterable_simple = lambda v: \
                False if hasattr(v, 'base_cls') else hasattr(v, '__iter__')

            getter = lambda v: \
                [None if int(w[1]) < 0 else self.stores[int(w[0])].load(int(w[1])) for w in v.tolist()] \
                    if len(v.shape) > 1 else None if int(v[1]) < 0 else self.stores[int(v[0])].load(int(v[1]))

            setter = lambda v: \
                np.array([(-1, -1) if w is None else self.save(w)[1:] for w in v], dtype=np.int32) \
                    if set_iterable_simple(v) else (-1, -1) if v is None else self.save(v)[1:]

        elif var_type.startswith('lazyobj.'):
            getter = lambda v: \
                [None if int(w) < 0 else LoaderProxy(store, int(w)) for w in v.tolist()] \
                    if get_is_iterable(v) else None if int(v) < 0 else LoaderProxy(store, int(v))
            setter = lambda v: \
                np.array([-1 if w is None else store.save(w) for w in v], dtype=np.int32) \
                    if set_is_iterable(v) else -1 if v is None else store.save(v)

        elif var_type == 'lazyobj':
            # arbitrary object

            set_iterable_simple = lambda v: \
                False if hasattr(v, 'base_cls') else hasattr(v, '__iter__')

            getter = lambda v: \
                [None if int(w[1]) < 0 else LoaderProxy(self.stores[int(w[0])],int(w[1])) for w in v.tolist()] \
                    if len(v.shape) > 1 else None if int(v[1]) < 0 else LoaderProxy(self.stores[int(v[0])],int(v[1]))

            setter = lambda v: \
                np.array([(-1, -1) if w is None else self.save(w)[1:] for w in v], dtype=np.int32) \
                    if set_iterable_simple(v) else (-1, -1) if v is None else self.save(v)[1:]

        elif var_type == 'store':
            setter = lambda v: v.prefix
            getter = lambda v: self.objects[v]

        return getter, setter, store

    def create_variable_delegate(self, var_name):
        """
        Create a delegate property that wraps the netcdf.Variable and takes care
        of type conversions
        """

        if var_name not in self.vars:
            var = self.variables[var_name]

            if not hasattr(var, 'var_type'):
                return

            getter, setter, store = self.create_type_delegate(var.var_type)

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
                        return lambda v: \
                            [None if hasattr(w, 'mask') else my_getter(w) for w in v] \
                                if type(v) is not str and len(v.shape) > 0 else \
                                (None if hasattr(v, 'mask') else my_getter(v))

                    if getter is not None:
                        getter = _get2(getter)
                    else:
                        getter = _get2(lambda v: v)

            self.vars[var_name] = NetCDFPlus.ValueDelegate(var, getter, setter, store)

        else:
            raise ValueError("Variable '%s' is already taken!")

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
            # last dimension is simply [] so we allow arbitrary length and remove the last dimension
            variable_length = True
            dimensions = dimensions[:-1]
        else:
            variable_length = False

        if var_type == 'obj' or var_type == 'lazyobj':
            dimensions.append('pair')
            if chunksizes is not None:
                chunksizes = tuple(list(chunksizes) + [2])

        nc_type = NetCDFPlus.var_type_to_nc_type(var_type)

        for dim_name, size in new_dimensions.items():
            ncfile.create_dimension(dim_name, size)

        dimensions = tuple(dimensions)

        # if chunksizes are strings then replace by the actual size of the dimension
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
            ncvar = ncfile.createVariable(var_name, vlen_t, dimensions,
                                          zlib=False, chunksizes=chunksizes)
        else:
            ncvar = ncfile.createVariable(var_name, nc_type, dimensions,
                                          zlib=False, chunksizes=chunksizes)

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
                raise NotImplementedError('Unit by abbreviated string representation is not yet supported')

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
                dim_names = map(lambda p: '#ix{0}:{1}'.format(*p), enumerate(dimensions))

            idx_desc = '[' + ']['.join(dim_names) + ']'
            description = var_name + idx_desc + ' is ' + description.format(idx=dim_names[0],
                                                                            ix=dim_names)

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

        if type(test_type) is np.ndarray:
            dimensions = test_type.shape
        else:
            if hasattr(test_value, '__len__'):
                dimensions = len(test_value)
                test_type = test_value[0]
                if NetCDFPlus.support_simtk_unit and type(test_type) is u.Quantity:
                    for val in test_value:
                        if type(val._value) is not type(test_value._value):
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
