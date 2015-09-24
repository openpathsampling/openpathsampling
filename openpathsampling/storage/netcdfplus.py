__author__ = 'jan-hendrikprinz'

import logging

logger = logging.getLogger(__name__)
init_log = logging.getLogger('openpathsampling.initialization')

from todict import StorableObjectJSON

from objproxy import LoaderProxy

import numpy as np
import netCDF4 as netcdf
import os.path

#=============================================================================================
# Extended NetCDF Storage for multiple forked trajectories
#=============================================================================================

class NetCDFPlus(netcdf.Dataset):
    support_simtk_unit = True

    class Value_Delegate(object):
        def __init__(self, variable, getter = None, setter = None):
            self.variable = variable

            if setter is None:
                setter = lambda v : v
            self.setter = setter

            if getter is None:
                getter = lambda v : v
            self.getter = getter

        def __setitem__(self, key, value):
            self.variable[key] = self.setter(value)

        def __getitem__(self, key):
            return self.getter(self.variable[key])

        def __getattr__(self, item):
            return getattr(self.variable, item)

    class Key_Delegate(object):
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
        return self._objects

    def _register_storages(self):
        pass

    def __init__(self, filename, mode=None, units=None):
        '''
        Create a storage for complex objects in a netCDF file

        Parameters
        ----------
        filename : string
            filename of the netcdf file to be used or created
        mode : string, default: None
            the mode of file creation, one of 'w' (write), 'a' (append) or
            None, which will append any existing files.
        template : openpathsampling.Snapshot
            a Snapshot instance that contains a reference to a Topology, the
            number of atoms and used units
        units : dict of {str : simtk.unit.Unit } or None
            representing a dict of string representing a dimension
            ('length', 'velocity', 'energy') pointing to
            the simtk.unit.Unit to be used. If not None overrides the
            standard units used
        '''

        if mode == None:
            if os.path.isfile(filename):
                logger.info("Open existing netCDF file '%s' for storage", filename)
                mode = 'a'
            else:
                logger.info("Create new netCDF file '%s' for storage", filename)
                mode = 'w'

        self.filename = filename

        # call netCDF4-python to create or open .nc file
        super(NetCDFPlus, self).__init__(filename, mode)

        self._setup_class()

        if units is not None:
            self.dimension_units.update(units)

        self._register_storages()
        self.set_caching_mode('default')

        if mode == 'w':
            logger.info("Setup netCDF file and create variables")

            # add shared dimension for everyone. scalar and spatial
            if 'scalar' not in self.dimensions:
                self.createDimension('scalar', 1) # scalar dimension

            self._initialize()

            logger.info("Finished setting up netCDF file")

        elif mode == 'a' or mode == 'r+' or mode == 'r':
            logger.debug("Restore the dict of units from the storage")

            # Create a dict of simtk.Unit() instances for all netCDF.Variable()
            for variable_name in self.variables:
                variable = self.variables[variable_name]

                if self.support_simtk_unit:
                    if hasattr(variable, 'unit_simtk'):
                        unit_dict = self.simplifier.from_json(getattr(variable, 'unit_simtk'))
                        if unit_dict is not None:
                            unit = self.simplifier.unit_from_dict(unit_dict)

                        self.units[str(variable_name)] = unit

            self.update_delegates()

            # After we have restored the units we can load objects from the storage
            self._restore()

        self.sync()

    def _setup_class(self):
        """
        Sets the basic properties for the storage
        """
        self._objects = {}
        self._storages = {}
        self._storages_base_cls = {}
        self._obj_store = {}
        self.simplifier = StorableObjectJSON(self)
        self.vars = dict()
        self.units = dict()

        self.dimension_units = dict()

    def add(self, name, store, register_attr=True):
        store.register(self, name)

        if register_attr:
            if hasattr(self, name):
                raise ValueError('Attribute name %s is already in use!' % name)

            setattr(self, store.prefix, store)


        self._objects[name] = store

        self._storages[store.content_class] = store
#        self._storages[store.content_class.__name__] = store
#        self._storages[store.content_class.__name__.lower()] = store

        self._obj_store[store.content_class] = store
        self._obj_store.update({cls : store for cls in store.content_class.descendants()})

    def _register(self):
        pass

    def _initialize(self):
        """
        Hook after a new file is created.

        This is used to setup all variables in the storage
        """
        pass

    def _restore(self):
        """
        Hook after an existing file is openend
        """
        pass

    def __repr__(self):
        return "Storage @ '" + self.filename + "'"

    def __getattr__(self, item):
        return self.__dict__[item]

    def __setattr__(self, key, value):
        self.__dict__[key] = value

    def _init_storages(self):
        '''
        Run the initialization on all added classes, when the storage is
        created only!

        Notes
        -----
        Only runs when the storage is created.
        '''

        for storage in self._objects.values():
            storage._init()

        self.update_delegates()

    def _initialize_netCDF(self):
        """
        Initialize the netCDF file for storage itself.
        """
        pass

    def list_stores(self):
        """
        Return a list of registered stores

        Returns
        -------
        list of str
            list of stores that can be accessed using `storage.[store]`
        """
        return [store.prefix for store in self.objects.values()]

    def list_storable_objects(self):
        """
        Return a list of storable object base classes

        Returns
        -------
        list of class
            list of base classes that can be stored using `storage.save(obj)`
        """
        return [store.content_class for store in self.objects.values()]

    def save(self, obj, *args, **kwargs):
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
            return [ self.save(part, *args, **kwargs) for part in obj]

        elif type(obj) is tuple:
            # a tuple will store all parts
            return tuple([self.save(part, *args, **kwargs) for part in obj])

        elif obj.__class__ in self._obj_store:
            # to store we just check if the base_class is present in the storages
            # also we assume that if a class has no base_cls
            store = self._obj_store[obj.__class__]
            store.save(obj, *args, **kwargs)

            # save has been called, all is good
            return True

        # Could not save this object. Might raise an exception, but
        # return an empty string as type
        raise RuntimeWarning("Objects of type '%s' cannot be stored!" %
                             obj.__class__.__name__)

    def load(self, obj_type, *args, **kwargs):
        """
        Load an object of the specified type from the storage

        Parameters
        ----------
        obj_type : str or class
            the store or class of the base object to be loaded.

        Returns
        -------
        object
            the object loaded from the storage

        Notes
        -----
        If you want to load a sub-classed Ensemble you need to load using
        `Ensemble` or `"Ensemble"` and not use the subclass
        """

        if obj_type in self._storages:
            store = self._storages[obj_type]
            return store.load(*args, **kwargs)
        elif obj_type in self._obj_store:
            # check if a store for the base_cls exists and use this one
            store = self._obj_store[obj_type]
            return store.load(*args, **kwargs)
        elif obj_type in self.simplifier.class_list:
            store = self._obj_store[self.simplifier.class_list[obj_type]]
            return store.load(*args, **kwargs)

        raise RuntimeError('No store registered to load variable type %s' % obj_type)

    def idx(self, obj):
        """
        Return the index used to store the given object in this storage

        Parameters
        ----------
        obj : object
            The stored object from which the index is to be returned
        """
        if hasattr(obj, 'base_cls'):
            store = self._storages[obj.base_cls]
            return store.idx(obj)

    def clone_storage(self, storage_to_copy, new_storage):
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
        if type(storage_to_copy) is str:
            storage_name = storage_to_copy
        else:
            storage_name = storage_to_copy.prefix

        copied_storages = 0

        for variable in self.variables.keys():
            if variable.startswith(storage_name + '_'):
                copied_storages += 1
                if variable not in new_storage.variables:
                    # collectivevariables have additional variables in the storage that need to be copied
                    # TODO: copy chunksizes?
                    new_storage.createVariable(variable, str(self.variables[variable].dtype), self.variables[variable].dimensions)
                    for attr in self.variables[variable].ncattrs():
                        setattr(new_storage.variables[variable], attr, getattr(self.variables[variable], attr))

                    new_storage.variables[variable][:] = self.variables[variable][:]
                else:
                    for idx in range(0, len(self.variables[variable])):
                        new_storage.variables[variable][idx] = self.variables[variable][idx]

        if copied_storages == 0:
            raise RuntimeWarning('Potential error in storage name. No storage variables copied from ' + storage_name)

    def create_dimension(self, dim_name, size = None):
        """
        Initialize a new dimension in the storage.
        Wraps the netCDF createDimension

        Parameters
        ----------
        name : str
            the name for the new dimension
        length : int
            the number of elements in this dimension. None (default) means
            an infinite dimension that extends when more objects are stored

        """
        if dim_name not in self.storage.dimensions:
            self.storage.createDimension(dim_name, size)

    def var_type_to_nc_type(self, var_type):
        type_conversion = {
            'float' : np.float32,
            'int' : np.int32,
            'long' : np.int64,
            'index' : np.int32,
            'length' : np.int32,
            'bool' : np.int16,
            'str' : 'str',
            'json' : 'str',
            'numpy.float32' : np.float32,
            'numpy.float64' : np.float32,
            'numpy.int8' : np.int8,
            'numpy.int16' : np.int16,
            'numpy.int32' : np.int32,
            'numpy.int64' : np.int64,
            'numpy.uint8' : np.uint8,
            'numpy.uinf16' : np.uint16,
            'numpy.uint32' : np.uint32,
            'numpy.uint64' : np.uint64,
            'obj' : np.int32
        }

        if var_type.startswith('obj.') or var_type.startswith('lazyobj.'):
            nc_type = np.int32
        else:
            nc_type = type_conversion[var_type]

        return nc_type

    def create_type_delegate(self, var_type):
        getter = None
        setter = None

        if var_type == 'int':
            getter = lambda v : v.tolist()
            setter = lambda v : np.array(v)

        elif var_type == 'bool':
            getter = lambda v : v.astype(np.bool).tolist()
            setter = lambda v : np.array(v, dtype=np.int8)

        elif var_type == 'index':
            getter = lambda v : \
                [ None if int(w) < 0 else int(w) for w in v.tolist()] \
                    if hasattr(v, '__iter__') else None if int(v) < 0 else int(v)
            setter = lambda v : \
                [ -1 if w is None else w for w in v] \
                    if hasattr(v, '__iter__') else -1 if v is None else v

        elif var_type == 'float':
            getter = lambda v : v.tolist()
            setter = lambda v : np.array(v)

        elif var_type.startswith('numpy.'):
            pass

        elif var_type == 'json':
            setter = lambda v : \
                self.simplifier.to_json(v)
            getter = lambda v : self.simplifier.from_json(v)

        elif var_type.startswith('obj.'):
            store = getattr(self, var_type[4:])
            base_type = store.content_class

            iterable = lambda v : \
                not v.base_cls is base_type if hasattr(v, 'base_cls') else hasattr(v, '__iter__')

            getter = lambda v : \
                [ None if int(w) < 0 else store.load(int(w)) for w in v.tolist()] \
                    if iterable(v) else None if int(v) < 0 else store.load(int(v))
            setter = lambda v : \
                np.array([ -1 if w is None else store.save(w) for w in v], dtype=np.int32) \
                    if iterable(v) else -1 if v is None else store.save(v)

        elif var_type.startswith('lazyobj.'):
            store = getattr(self, var_type[8:])
            base_type = store.content_class

            iterable = lambda v : \
                not v.base_cls is base_type if hasattr(v, 'base_cls') else hasattr(v, '__iter__')

            getter = lambda v : \
                [ None if int(w) < 0 else LoaderProxy({store : int(w)}) for w in v.tolist()] \
                    if iterable(v) else None if int(v) < 0 else LoaderProxy({store : int(v)})
            setter = lambda v : \
                np.array([ -1 if w is None else store.save(w) for w in v], dtype=np.int32) \
                    if iterable(v) else -1 if v is None else store.save(v)

        return getter, setter

    def create_variable_delegate(self, var_name):
        """
        Create a delegate property that wraps the netcdf.Variable and takes care
        of type conversions
        """

        if var_name not in self.vars:
            var = self.variables[var_name]

            if not hasattr(var, 'var_type'):
                return

            getter, setter = self.create_type_delegate(var.var_type)

            if True or self.support_simtk_unit:
                if hasattr(var, 'unit_simtk'):
                    if var_name not in self.units:
                        self.update_simtk_unit(var_name)

                    unit = self.units[var_name]

                    def _get(my_getter):
                        import simtk.unit as u
                        if my_getter is None:
                            return  lambda v : u.Quantity(v, unit)
                        else:
                            return  lambda v : u.Quantity(my_getter(v), unit)
                    def _set(my_setter):
                        if my_setter is None:
                            return  lambda v : v / unit
                        else:
                            return  lambda v : my_setter(v / unit)
                    getter = _get(getter)
                    setter = _set(setter)

            if True:
                if hasattr(var, 'maskable'):
                    iterable = lambda v : type(v) is not int and hasattr(v, '__iter__')

                    def _get2(my_getter):
                        return lambda v : \
                            [None if hasattr(w, 'mask') else my_getter(w) for w in v ] \
                            if type(v) is not str and len(v.shape) > 0 else \
                                (None if hasattr(v, 'mask') else my_getter(v))

                    if getter is not None:
                        getter = _get2(getter)
                    else:
                        getter = _get2(lambda v : v)

            self.vars[var_name] = NetCDFPlus.Value_Delegate(var, getter, setter)

        else:
            raise ValueError("Variable '%s' is already taken!")

    def create_variable(self, var_name,
                        var_type,
                        dimensions,
                        description=None,
                        chunksizes=None,
                        simtk_unit=None,
                        maskable=False,
                        index_by=None,
    ):
        '''
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
        '''

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

        nc_type = self.var_type_to_nc_type(var_type)

        for dim_name, size in new_dimensions.items():
            ncfile.create_dimension(dim_name, size)

        dimensions = tuple(dimensions)

        if variable_length:
            vlen_t = ncfile.createVLType(nc_type, var_name + '_vlen')
            ncvar = ncfile.createVariable(var_name, vlen_t, dimensions,
                                          zlib=False, chunksizes=chunksizes)
        else:
            ncvar = ncfile.createVariable(var_name, nc_type, dimensions,
                                          zlib=False, chunksizes=chunksizes)

        setattr(ncvar,      'var_type', var_type)

        if self.support_simtk_unit and simtk_unit is not None:

            import simtk.unit as u

            unit_instance = u.Unit({})
            symbol = 'none'

            if isinstance(simtk_unit, u.Unit):
                unit_instance = simtk_unit
                symbol = unit_instance.get_symbol()
            elif isinstance(simtk_unit, u.BaseUnit):
                unit_instance = u.Unit({simtk_unit : 1.0})
                symbol = unit_instance.get_symbol()
            elif type(simtk_unit) is str and hasattr(u, simtk_unit):
                unit_instance = getattr(u, simtk_unit)
                symbol = unit_instance.get_symbol()
            elif type(simtk_unit) is str and simtk_unit in self.dimension_units:
                unit_instance = self.dimension_units[simtk_unit]
                symbol = unit_instance.get_symbol()
            else:
                raise NotImplementedError('Unit by abbreviated string representation is not yet supported')

            json_unit = self.simplifier.unit_to_json(unit_instance)


            # store the unit in the dict inside the Storage object
            self.units[var_name] = unit_instance

            # Define units for a float variable
            setattr(ncvar,      'unit_simtk', json_unit)
            setattr(ncvar,      'unit', symbol)

        if maskable:
            setattr(ncvar,      'maskable', 'True')

        if description is not None:
            if type(dimensions) is str:
                dim_names = [ dimensions ]
            else:
                dim_names = map( lambda p : '#ix{0}:{1}'.format(*p), enumerate(dimensions))

            idx_desc = '[' + ']['.join(dim_names) + ']'
            description = var_name + idx_desc + ' is ' + description.format(idx=dim_names[0], ix=dim_names)

            # Define long (human-readable) names for variables.
            setattr(ncvar,    "long_str", description)

        return ncvar

    def update_delegates(self):
        for name in self.variables:
            if name not in self.vars:
                self.create_variable_delegate(name)
