'''
Created on 06.07.2014

@author: JDC Chodera
@author: JH Prinz
'''

import netCDF4 as netcdf
import os.path

import logging

logger = logging.getLogger(__name__)
init_log = logging.getLogger('openpathsampling.initialization')

import numpy as np
import simtk.unit as u
import re
import functools

import openpathsampling as paths
from openpathsampling.todict import ObjectJSON

from openpathsampling.tools import LRUCache, WeakLRUCache, WeakCache, WeakLimitCache

#=============================================================================================
# SOURCE CONTROL
#=============================================================================================

__version__ = "$Id: NoName.py 1 2014-07-06 07:47:29Z jprinz $"

#=============================================================================================
# NetCDF Storage for multiple forked trajectories
#=============================================================================================

class NetCDFPlus(netcdf.Dataset):
    '''
    A netCDF4 wrapper to store trajectories based on snapshots of an OpenMM
    simulation. This allows effective storage of shooting trajectories '''

    # I was wondering if it makes sense to just treat the stores as a
    # variable with one index, since their basic behaviour is the same.
    # For now I would leave it this way to make clear, that means a normal
    # netCDF variable is accessed by storage.variables[name][idx] and the
    # `special` variable is accessed by storage.store[idx]

    allowed_types = [
        'int', 'float', 'long', 'str', 'bool',
        'numpy.float32', 'numpy.float64',
        'numpy.int8', 'numpy.inf16', 'numpy.int32', 'numpy.int64',
        'numpy.uint8', 'numpy.uinf16', 'numpy.uint32', 'numpy.uint64',
        'index', 'json'
    ]
    # Default caching for online use


    @staticmethod
    def default_cache_sizes():
        return {
            'trajectories' : WeakLRUCache(10000),
            'snapshots' : WeakLRUCache(50000),
            'configurations' : WeakLRUCache(10000),
            'momentum' : WeakLRUCache(10000),
            'samples' : WeakLRUCache(25000),
            'samplesets' : False,
            'collectivevariables' : True,
            'pathmovers' : True,
            'shootingpoints' : WeakLRUCache(10000),
            'shootingpointselectors' : True,
            'engines' : True,
            'pathsimulators' : True,
            'volumes' : True,
            'ensembles' : True,
            'pathmovechanges' : False,
            'transitions' : True,
            'networks' : True,
            '_details' : False,
            'steps' : WeakLRUCache(1000)
        }

    # Analysis caching is very large to allow fast processing

    @staticmethod
    def analysis_cache_sizes():
        return {
            'trajectories' : WeakLimitCache(100000),
            'snapshots' : WeakLimitCache(500000),
            'configurations' : WeakLRUCache(10000),
            'momentum' : WeakLRUCache(10000),
            'samples' : WeakLimitCache(250000),
            'samplesets' : WeakLimitCache(100000),
            'collectivevariables' : True,
            'pathmovers' : True,
            'shootingpoints' : WeakLimitCache(100000),
            'shootingpointselectors' : True,
            'engines' : True,
            'pathsimulators' : True,
            'volumes' : True,
            'ensembles' : True,
            'pathmovechanges' : WeakLimitCache(250000),
            'transitions' : True,
            'networks' : True,
            '_details' : False,
            'steps' : WeakLimitCache(100000)
        }

    # Production. No loading, only last 1000 steps and a few other objects for error
    # testing

    def production_cache_sizes(self):
        return {
            'trajectories' : WeakLRUCache(),
            'snapshots' : WeakLRUCache(),
            'configurations' : WeakLRUCache(),
            'momentum' : WeakLRUCache(),
            'samples' : WeakLRUCache(),
            'samplesets' : False,
            'collectivevariables' : False,
            'pathmovers' : False,
            'shootingpoints' : False,
            'shootingpointselectors' : False,
            'engines' : False,
            'pathsimulators' : False,
            'volumes' : False,
            'ensembles' : False,
            'pathmovechanges' : False,
            'transitions' : False,
            'networks' : False,
            '_details' : False,
            'steps' : WeakCache()
        }

    # No caching (so far only CVs internal storage is there)

    def no_cache_sizes(self):
        return {
            'trajectories' : False,
            'snapshots' : False,
            'configurations' : False,
            'momentum' : False,
            'samples' : False,
            'samplesets' : False,
            'collectivevariables' : False,
            'pathmovers' : False,
            'shootingpoints' : False,
            'shootingpointselectors' : False,
            'engines' : False,
            'pathsimulators' : False,
            'volumes' : False,
            'ensembles' : False,
            'pathmovechanges' : False,
            'transitions' : False,
            'networks' : False,
            '_details' : False,
            'steps' : False
        }

    class Variable_Delegate(object):
        def __init__(self, variable, getter, setter):
            self.variable = variable

            if setter is None:
                setter = lambda v : v
            self.setter = setter

            if getter is None:
                getter = lambda v : v
            self.getter = getter

        def __setitem__(self, key, value):
            self.variable[key] = self.setter(value)

        def __getitem__(self, item):
            return self.getter(self.variable[item])

        def __getattr__(self, item):
            return getattr(self.variable, item)

    def set_caching_mode(self, mode='default'):
        """
        Set default values for all caches

        Parameters
        ----------
        caching : str
            One of the following values is allowed 'default', 'production',
            'analysis', 'off'

        """

        available_cache_sizes = {
            'default': self.default_cache_sizes,
            'analysis': self.analysis_cache_sizes,
            'production': self.production_cache_sizes,
            'off': self.no_cache_sizes
        }

        if mode in available_cache_sizes:
            # We need cache sizes as a function. Otherwise we will reuse the same
            # caches for each storage and that will cause problems! Lots of...
            cache_sizes = available_cache_sizes[mode]()
        else:
            raise ValueError(
                "mode '" + mode + "' is not supported. Try one of " +
                str(available_cache_sizes.keys())
            )

        for store_name, caching in cache_sizes.iteritems():
            if hasattr(self, store_name):
                store = getattr(self, store_name)
                store.set_caching(caching)

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
                unit = None
                variable = self.variables[variable_name]

                if hasattr(variable, 'unit_simtk'):
                    unit_dict = self.simplifier.from_json(getattr(variable, 'unit_simtk'))
                    if unit_dict is not None:
                        unit = self.simplifier.unit_from_dict(unit_dict)

                self.units[str(variable_name)] = unit

            self.update_delegates()

            # After we have restored the units we can load objects from the storage
            self._restore()

        self.sync()

        self.set_caching_mode('default')

    @property
    def objects(self):
        return self._objects

    def _setup_class(self):
        """
        Sets the basic properties for the storage
        """
        self._objects = {}
        self._storages = {}
        self._storages_base_cls = {}
        self._obj_store = {}
        self.simplifier = StorableObjectJSON(self)
        self.units = dict()
        self.vars = dict()
        # use no units
        self.dimension_units = {
            'length': u.Unit({}),
            'velocity': u.Unit({}),
            'energy': u.Unit({})
        }
        # use MD units
        self.dimension_units = {
            'length': u.nanometers,
            'velocity': u.nanometers / u.picoseconds,
            'energy': u.kilojoules_per_mole
        }

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

    def get_unit(self, dimension):
        """
        Return a simtk.Unit instance from the unit_system the is of the specified dimension, e.g. length, time
        """
        return u.Unit({self.unit_system.base_units[u.BaseDimension(dimension)] : 1.0})

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
        # no store found. This is bad and should be logged and raise
        # an exception
#        return None

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

    def create_dimension(self, dimname, size = None):
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
        if dimname not in self.storage.dimensions:
            self.storage.createDimension(dimname, size)

    @staticmethod
    def compose(*functions):
        funcs = [func for func in functions if func is not None]
        return functools.reduce(lambda f, g: lambda x: f(g(x)), funcs, lambda x: x)


    def parse_var_type(self, var_name, var_type):
        dimensions = None
        new_dimensions = dict()
        if '[' in var_type:
            # find dimensions
            dimensions = re.findall('\[([a-zA-Z]*|[0-9]*)\]', var_type)

        for ix, dim in enumerate(dimensions):
            try:
                length = int(dim)
                dimensions[ix] = var_name + '_dim_' + str(ix)
                new_dimensions[dimensions[ix]] = length
            except ValueError:
                pass

        if dimensions[-1] == '':
            # last dimension is simply [] so we allow arbitrary length and remove the last dimension
            variable_length = True
            dimensions = dimensions[:-1]
        else:
            variable_length = False

        inner_type = None
        outer_fncs = []

        if '(' in var_type:
            red_var_type = var_type
            first_group = True
            while '(' in red_var_type:
                outer_params = list(re.search('([a-zA-Z\._]+)\([ ]*(?:([a-zA-Z0-9_\.]+)[ ]*,)*[ ]*([a-zA-Z0-9_\.]+)[ ]*\)', red_var_type).groups())
                outer_type = outer_params[0]
                if first_group:
                    if outer_params[1] is None:
                        inner_type = outer_params[2]
                        outer_params = [outer_params[0]] + outer_params[2:]
                    else:
                        inner_type = outer_params[1]
                    first_group = False
                outer_args = list(outer_params)[2:]
                outer_fncs.append([outer_type, outer_args])
                red_var_type = re.sub('(([a-zA-Z\._]+)\([ ]*(?:([a-zA-Z0-9_\.]+)[ ]*,)*[ ]*([a-zA-Z0-9_\.]+)[ ]*\))', 'inner', red_var_type)
        else:
            inner_type = re.search('[ ]*([a-zA-Z0-9_\.]+)[ ]*', var_type).group()

        inner = self.create_type_delegate(inner_type)

        if not inner['atomic']:
            raise ValueError('Type of inner function must be numeric!')

        get_list = [inner['getter']]
        set_list = [inner['setter']]

        for outer_type, outer_args in outer_fncs:
            outer = self.create_type_delegate(outer_type, *outer_args)

            if outer['atomic']:
                raise ValueError('Type of outer function must be')
            _get = outer['getter']
            _set = outer['setter']

            get_list.insert(0,_get)
            set_list.append(_set)

        getter = self.compose(*get_list)
        setter = self.compose(*set_list)

        type_conversion = {
            'float' : np.float32,
            'int' : np.int32,
            'long' : np.int64,
            'index' : np.int32,
            'length' : np.int32,
            'bool' : np.uint8,
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
            'numpy.uint64' : np.uint64
        }

        if inner_type.startswith('obj.'):
            nc_type = np.int32
        else:
            nc_type = type_conversion[inner_type]

        return {
            'getter' : getter,
            'setter' : setter,
            'dimensions' : tuple(dimensions),
            'new_dimensions' : new_dimensions,
            'variable_length' : variable_length,
            'nc_type' : nc_type,
            'outer' : outer_fncs,
            'inner' : inner_type,
            'get_list' : get_list,
            'set_list' : set_list
        }


    def create_type_delegate(self, var_type, *args):
        getter = None
        setter = None
        n_args = 0
        atomic = True

        if var_type == 'int':
            getter = lambda v : v.tolist()
            setter = lambda v : np.array(v)

        elif var_type == 'bool':
            getter = lambda v : v.astype(np.bool).tolist()
            setter = lambda v : np.array(v, dtype=np.uint8)

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

        elif var_type == 'tuple':
            atomic = False
            getter = tuple

        elif var_type == 'list':
            atomic = False
            getter = list # might be obsolete

        elif var_type == 'json':
            atomic = True
            setter = lambda v : \
                self.simplifier.to_json(v)
            getter = lambda v : self.simplifier.from_json(v)

        elif var_type == 'quantity':
            atomic = False
            n_args = 1
            unit = self.dimension_units[args[0]]
            getter = lambda v : u.Quantity(v, unit)
            setter = lambda v : v / unit

        elif var_type.startswith('obj.'):
            atomic = True
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

        return {
            'getter' : getter,
            'setter' : setter,
            'atomic' : atomic,
            'n_args' : n_args
        }

    def create_variable_delegate(self, name):
        """
        Create a delegate property that wraps the netcdf.Variable and takes care
        of type conversions
        """

        var = self.variables[name]

        if not hasattr(var, 'var_type'):
            return

        var_info = self.parse_var_type(name, var.var_type)

        getter = var_info['getter']
        setter = var_info['setter']

        self.vars[name] = NetCDFPlus.Variable_Delegate(var, getter, setter)

    def create_variable(self, name, var_type,
                      description=None, chunksizes=None):
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
        var_info = self.parse_var_type(name, var_type)

        dimensions = var_info['dimensions']
        variable_length = var_info['variable_length']
        nc_type = var_info['nc_type']

        for dimname, size in var_info['new_dimensions'].items():
            ncfile.create_dimension(dimname, size)

        if variable_length:
            vlen_t = ncfile.createVLType(nc_type, name + '_vlen')
            ncvar = ncfile.createVariable(name, vlen_t, dimensions,
                                          zlib=False, chunksizes=chunksizes)
        else:
            ncvar = ncfile.createVariable(name, nc_type, dimensions,
                                          zlib=False, chunksizes=chunksizes)

        setattr(ncvar,      'var_type', var_type)

        # if var_type == 'float' or units is not None:
        #
        #     unit_instance = u.Unit({})
        #     symbol = 'none'
        #
        #     if isinstance(units, u.Unit):
        #         unit_instance = units
        #         symbol = unit_instance.get_symbol()
        #     elif isinstance(units, u.BaseUnit):
        #         unit_instance = u.Unit({units : 1.0})
        #         symbol = unit_instance.get_symbol()
        #     elif type(units) is str and hasattr(u, units):
        #         unit_instance = getattr(u, units)
        #         symbol = unit_instance.get_symbol()
        #     elif type(units) is str and units is not None:
        #         symbol = units

#            json_unit = self.simplifier.unit_to_json(unit_instance)

#            # store the unit in the dict inside the Storage object
#            self.units[name] = unit_instance
#
#            # Define units for a float variable
#            setattr(ncvar,      'unit_simtk', json_unit)
#            setattr(ncvar,      'unit', symbol)

        if description is not None:
            if type(dimensions) is str:
                dim_names = [ dimensions ]
            else:
                dim_names = map( lambda p : '#ix{0}:{1}'.format(*p), enumerate(dimensions))

            idx_desc = '[' + ']['.join(dim_names) + ']'
            description = name + idx_desc + ' is ' + description.format(idx=dim_names[0], ix=dim_names)

            # Define long (human-readable) names for variables.
            setattr(ncvar,    "long_str", description)

    def update_delegates(self):
        for name in self.variables:
            if name not in self.vars:
                self.create_variable_delegate(name)


class Storage(NetCDFPlus):

    @property
    def template(self):
        """
        Return the template snapshot from the storage

        Returns
        -------
        Snapshot
            the initial snapshot
        """
        if self._template is None:
            self._template = self.snapshots.load(int(self.variables['template_idx'][0]))

        return self._template

    def clone(self, filename, subset):
        """
        Creates a copy of the netCDF file and allows to reduce the used atoms.

        Notes
        -----
        This is mostly used to remove water but keep the data intact.
        """

        storage2 = Storage(filename=filename, template=self.template.subset(subset), mode='w')

        # Copy all configurations and momenta to new file in reduced form
        for obj in self.configurations:
            storage2.configurations.save(obj.copy(subset), idx=obj.idx[self.configurations])
        for obj in self.momenta.iterator():
            storage2.momenta.save(obj.copy(subset), idx=obj.idx[self.momenta])

        # All other should be copied one to one. We do this explicitely although we could just copy all
        # and exclude configurations and momenta, but this seems cleaner

        for storage_name in [
                'pathmovers', 'topologies', 'networks', 'details', 'trajectories',
                'shootingpoints', 'shootingpointselectors', 'engines', 'volumes',
                'samplesets', 'ensembles', 'transitions', 'steps', 'pathmovechanges',
                'samples', 'snapshots', 'pathsimulators', 'cvs'
            ]:
            self.clone_storage(storage_name, storage2)

        storage2.close()

    # TODO: Need to copy cvs without caches!
    def clone_empty(self, filename):
        """
        Creates a copy of the netCDF file and replicates only the static parts which I consider
            ensembles, volumes, engines, pathmovers, shootingpointselectors. We do not need to
            reconstruct collectivevariables since these need to be created again completely and then
            the necessary arrays in the file will be created automatically anyway.

        Notes
        -----
        This is mostly used to restart with a fresh file. Same setup, no results.
        """
        storage2 = Storage(filename=filename, template=self.template, mode='w')

        for storage_name in [
                'pathmovers', 'topologies', 'networks',
                'shootingpointselectors', 'engines', 'volumes',
                'ensembles', 'transitions', 'pathsimulators'
            ]:
            self.clone_storage(storage_name, storage2)

        storage2.close()

    @property
    def n_atoms(self):
        return self.topology.n_atoms

    @property
    def n_spatial(self):
        return self.topology.n_spatial

    def __init__(self, filename, mode=None,
                 template=None, units=None):
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

        self._template = template
        super(Storage, self).__init__(filename, mode, units=units)

    def _register_storages(self):
        """
        Register all Stores used in the OpenPathSampling Storage

        """

        # objects with special storages

        # self.objectname = ... could also be done in the initialization
        # automatically. But the IDE would not be able to autocomplete
        # so we leave it this way :)

        self.add('trajectories', paths.storage.TrajectoryStore())
        self.add('snapshots', paths.storage.SnapshotStore())
        self.add('configurations', paths.storage.ConfigurationStore())
        self.add('momenta', paths.storage.MomentumStore())
        self.add('samples', paths.storage.SampleStore())
        self.add('samplesets', paths.storage.SampleSetStore())
        self.add('pathmovechanges', paths.storage.PathMoveChangeStore())
        self.add('steps', paths.storage.MCStepStore())

        self.add('cvs', paths.storage.ObjectDictStore(paths.CollectiveVariable, paths.Snapshot))
        self.collectivevariables = self.cvs

        # normal objects

        self.add('details', paths.storage.ObjectStore(paths.Details, has_uid=False, has_name=False))
        self.add('topologies', paths.storage.ObjectStore(paths.Topology, has_uid=True, has_name=True))
        self.add('pathmovers', paths.storage.ObjectStore(paths.PathMover, has_uid=True, has_name=True))
        self.add('shootingpoints', paths.storage.ObjectStore(paths.ShootingPoint, has_uid=False))
        self.add('shootingpointselectors', paths.storage.ObjectStore(paths.ShootingPointSelector, has_uid=False, has_name=True))
        self.add('engines', paths.storage.ObjectStore(paths.DynamicsEngine, has_uid=True, has_name=True))
        self.add('pathsimulators', paths.storage.ObjectStore(paths.PathSimulator, has_uid=True, has_name=True))
        self.add('transitions', paths.storage.ObjectStore(paths.Transition, has_uid=True, has_name=True))
        self.add('networks', paths.storage.ObjectStore(paths.TransitionNetwork, has_uid=True, has_name=True))

        # nestable objects

        self.add('volumes', paths.storage.ObjectStore(paths.Volume, has_uid=True, nestable=True, has_name=True))
        self.add('ensembles', paths.storage.ObjectStore(paths.Ensemble, has_uid=True, nestable=True, has_name=True))

    def _initialize(self):
        # Set global attributes.
        setattr(self, 'title', 'OpenPathSampling Storage')
        setattr(self, 'ConventionVersion', '0.2')

        template = self._template

        if template.topology is not None:
            self.topology = template.topology
        else:
            raise RuntimeError("A Storage needs a template snapshot with a topology")

        if 'atom' not in self.dimensions:
            self.createDimension('atom', self.topology.n_atoms)

        # spatial dimensions
        if 'spatial' not in self.dimensions:
            self.createDimension('spatial', self.n_spatial)

        # update the units for dimensions from the template
        self.dimension_units.update(paths.tools.units_from_snapshot(template))

        self._init_storages()

        logger.info("Saving topology")

        topology_id = self.topologies.save(self.topology)

        logger.info("Create initial template snapshot")

        # Save the initial configuration
        self.snapshots.save(template)

        self.createVariable('template_idx', 'i4', 'scalar')
        self.variables['template_idx'][:] = template.idx[self.snapshots]

    def _restore(self):
        self.topology = self.topologies[0]


class StorableObjectJSON(ObjectJSON):
    def __init__(self, storage, unit_system = None):
        super(StorableObjectJSON, self).__init__(unit_system)
        self.excluded_keys = ['idx', 'json', 'identifier']
        self.storage = storage

    def simplify(self,obj, base_type = ''):
        if obj is self.storage:
            return { '_storage' : 'self' }
        if obj.__class__.__module__ != '__builtin__':
            if obj.__class__ in self.storage._obj_store:
                store = self.storage._obj_store[obj.__class__]
                if not store.nestable or obj.base_cls_name != base_type:
                    # this also returns the base class name used for storage
                    # store objects only if they are not creatable. If so they will only be created in their
                    # top instance and we use the simplify from the super class ObjectJSON
                    self.storage.save(obj)
                    return { '_idx' : obj.idx[store], '_obj' : store.prefix }

        return super(StorableObjectJSON, self).simplify(obj, base_type)

    def build(self,obj):
        if type(obj) is dict:
            if '_storage' in obj:
                if obj['_storage'] == 'self':
                    return self.storage

            if '_idx' in obj and '_obj' in obj:
                store = self.storage._objects[obj['_obj']]
                result = store.load(obj['_idx'])

                return result

        return super(StorableObjectJSON, self).build(obj)

class AnalysisStorage(Storage):
    def __init__(self, filename):
        super(AnalysisStorage, self).__init__(
            filename=filename,
            mode='r'
        )

        self.set_caching_mode('analysis')

        # Let's go caching

        self.samples.cache_all()
        self.samplesets.cache_all()
        self.cvs.cache_all()
        map(lambda x : x.cache_all(self.cvs), self.cvs)
        self.volumes.cache_all()
        self.ensembles.cache_all()
        self.pathmovers.cache_all()
        self.pathmovechanges.cache_all()
        self.steps.cache_all()
        