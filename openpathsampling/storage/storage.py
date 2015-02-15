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

import numpy
import simtk.unit as u

import openpathsampling as paths
import openpathsampling.todict

#=============================================================================================
# SOURCE CONTROL
#=============================================================================================

__version__ = "$Id: NoName.py 1 2014-07-06 07:47:29Z jprinz $"

#=============================================================================================
# NetCDF Storage for multiple forked trajectories
#=============================================================================================

class Storage(netcdf.Dataset):
    '''
    A netCDF4 wrapper to store trajectories based on snapshots of an OpenMM
    simulation. This allows effective storage of shooting trajectories '''

    # I was wondering if it makes sense to just treat the stores as a
    # variable with one index, since their basic behaviour is the same.
    # For now I would leave it this way to make clear, that means a normal
    # netCDF variable is accessed by storage.variables[name][idx] and the
    # `special` variable is accessed by storage.store[idx]

    def _register_storages(self, storage = None):
        """
        Register all Stores used in the OpenPathSampling Storage

        Parameters
        ----------
        storage : Storage
            the Storage object the store should use. Usually (default) the
            storage itself

        """
        if storage is None:
            storage = self

        # objects with special storages

        self.trajectory = paths.storage.TrajectoryStore(storage)
        self.snapshot = paths.storage.SnapshotStore(storage)
        self.configuration = paths.storage.ConfigurationStore(storage)
        self.momentum = paths.storage.MomentumStore(storage)
        self.sample = paths.storage.SampleStore(storage)
        self.sampleset = paths.storage.SampleSetStore(storage)

        self.collectivevariable = paths.storage.ObjectDictStore(storage, paths.OrderParameter, paths.Configuration)
        self.cv = self.collectivevariable

        # normal objects

        self.pathmover = paths.storage.ObjectStore(storage, paths.PathMover, is_named=True)
        self.movedetails = paths.storage.ObjectStore(storage, paths.MoveDetails, is_named=False)
        self.movepaths = paths.storage.ObjectStore(storage, paths.MovePath, is_named=False, nestable=True)
        self.shootingpoint = paths.storage.ObjectStore(storage, paths.ShootingPoint, is_named=False)
        self.shootingpointselector = paths.storage.ObjectStore(storage, paths.ShootingPointSelector, is_named=False)
        self.engine = paths.storage.ObjectStore(storage, paths.DynamicsEngine, is_named=True)

        # nestable objects

        self.volume = paths.storage.ObjectStore(storage, paths.Volume, is_named=True, nestable=True)
        self.ensemble = paths.storage.ObjectStore(storage, paths.Ensemble, is_named=True, nestable=True)

    def _setup_class(self):
        """
        Sets the basic properties for the storage
        """
        self._storages = {}
        self._storages_base_cls = {}
        self.links = []
        self.simplifier = paths.ObjectJSON()
        self.units = dict()
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

        if mode == None:
            if os.path.isfile(filename):
                logger.info("Open existing netCDF file '%s' for storage", filename)
                mode = 'a'
            else:
                logger.info("Create new netCDF file '%s' for storage", filename)
                mode = 'w'


        self.filename = filename

        super(Storage, self).__init__(filename, mode)

        self._setup_class()

        if units is not None:
            self.dimension_units.update(units)

        self._register_storages()

        if mode == 'w':
            logger.info("Setup netCDF file and create variables")

            if template.topology is not None:
                self.topology = template.topology
            else:
                raise RuntimeError("Storage need a template snapshot with topology")

            self._initialize_netCDF()

            # update the units for dimensions from the template
            self.dimension_units.update(paths.tools.units_from_snapshot(template))
            self._init_storages()

            logger.info("Saving topology")

            # create a json from the mdtraj.Topology() and store it
            self.write_str('topology', self.simplifier.to_json(self.topology))


            logger.info("Create initial template snapshot")

            # Save the initial configuration
            self.snapshot.save(template)

            self.createVariable('template_idx', 'i4', 'scalar')
            self.variables['template_idx'][:] = template.idx[self]

            self.sync()

            logger.info("Finished setting up netCDF file")

        elif mode == 'a' or mode == 'r+' or mode == 'r':
            self._restore_storages()

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

            # After we have restored the units we can load objects from the storage

            self.topology = self.simplifier.from_json(self.variables['topology'][0])

    def __repr__(self):
        return "OpenPathSampling netCDF Storage @ '" + self.filename + "'"

    @property
    def n_atoms(self):
        return self.topology.n_atoms

    @property
    def n_spatial(self):
        return self.topology.n_spatial

    @property
    def template(self):
        """
        Return the template snapshot from the storage

        Returns
        -------
        Snapshot
            the initial snapshot
        """
        return self.snapshot.load(int(self.variables['template_idx'][0]))

    def get_unit(self, dimension):
        """
        Return a simtk.Unit instance from the unit_system the is of the specified dimension, e.g. length, time
        """
        return u.Unit({self.unit_system.base_units[u.BaseDimension(dimension)] : 1.0})

    def __getattr__(self, item):
#        if item in self.__dict__:
            return self.__dict__[item]
 #       else:
  #          return super(Storage, self).__getattr__(item)

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

        for storage in self.links:
            # create a member variable which is the associated Class itself
            storage.dimension_units.update(units=self.dimension_units)
            storage._init()

    def _restore_storages(self):
        '''
        Run restore on all added classes. Usually there is nothing to do.
        '''
#        for storage in self.links:
#            storage._restore()
        pass

    def _initialize_netCDF(self):
        """
        Initialize the netCDF file for storage itself.
        """

        # add shared dimension for everyone. scalar and spatial
        if 'scalar' not in self.dimensions:
            self.createDimension('scalar', 1) # scalar dimension

        if 'atom' not in self.dimensions:
            self.createDimension('atom', self.topology.n_atoms)

        # spatial dimensions
        if 'spatial' not in self.dimensions:
            self.createDimension('spatial', self.n_spatial)

        # Set global attributes.
        setattr(self, 'title', 'Open-Transition-Interface-Sampling')
        setattr(self, 'application', 'Host-Guest-System')
        setattr(self, 'program', 'run.py')
        setattr(self, 'programVersion', __version__)
        setattr(self, 'Conventions', 'Multi-State Transition Interface TPS')
        setattr(self, 'ConventionVersion', '0.1')

        # Create a string to hold the topology
        self.init_str('topology')

        # Force sync to disk to avoid data loss.
        self.sync()

    def write_str(self, name, string):
        '''
        Write a string into a netCDF variable

        Parameters
        ----------
        name : str
            the name of the variable
        string : str
            the string to store
        '''
        packed_data = numpy.empty(1, 'O')
        packed_data[0] = string
        self.variables[name][:] = packed_data

    def load_str(self, name):
        '''
        Load a string from a netCDF variable

        Parameters
        ----------
        name : str
            the name of the variable

        Returns
        -------
        str
            the loaded str

        '''
        return self.variables[name][0]


    def init_str(self, name):
        '''
        Initialize a netCDF Variable to store a single string

        Parameters
        ----------
        name : str
            the name of the variable
        '''
        self.createVariable(name, 'str', 'scalar')

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
        elif hasattr(obj, 'base_cls'):
            # to store we just check if the base_class is present in the storages
            # also we assume that if a class has no base_cls
            if obj.base_cls in self._storages:
                store = self._storages[obj.base_cls]
                store.save(obj, *args, **kwargs)

                # save has been called, all is good
                return True

        # Could not save this object. Might raise an exception, but
        # return an empty string as type
        return False

    def load(self, obj_type, *args, **kwargs):
        """
        Load an object of the specified type from the storage

        Parameters
        ----------
        obj_type : str or class
            the string or class of the base object to be loaded.

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
        elif hasattr(obj_type, 'base_cls'):
            # check if a store for the base_cls exists and use this one
            if obj_type.base_cls in self._storages:
                store = self._storages[obj_type.base_cls]
                return store.load(*args, **kwargs)

        # no store found. This is bad and should be logged and raise
        # an exception
        return None

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


    def clone(self, filename, subset):
        """
        Creates a copy of the netCDF file and allows to reduce the used atoms.

        Notes
        -----
        This is mostly used to remove water but keep the data intact.
        """
        storage2 = Storage(filename=filename, template=self.template.subset(subset), mode='w')

        # Copy all configurations and momenta to new file in reduced form

        for obj in self.configuration.iterator():
            storage2.configuration.save(obj.copy(subset), idx=obj.idx[self])
        for obj in self.momentum.iterator():
            storage2.momentum.save(obj.copy(subset), idx=obj.idx[self])

        # All other should be copied one to one. We do this explicitely although we could just copy all
        # and exclude configurations and momenta, but this seems cleaner

        for storage_name in [
                'trajectory', 'snapshot', 'sample', 'sampleset', 'orderparameter',
                'pathmover', 'engine', 'movedetails', 'shootingpoint', 'shootingpointselector',
                'globalstate', 'volume', 'ensemble' ]:
            self.clone_storage(storage_name, storage2)

        storage2.close()

    def clone_empty(self, filename):
        """
        Creates a copy of the netCDF file and replicates only the static parts which I consider
            ensembles, volumes, engines, pathmovers, shootingpointselectors. We do not need to
            reconstruct orderparameters since these need to be created again completely and then
            the necessary arrays in the file will be created automatically anyway.

        Notes
        -----
        This is mostly used to restart with a fresh file. Same setup, no results.
        """
        storage2 = Storage(filename=filename, template=self.template, mode='w')

        for storage_name in [
                'pathmover', 'engine', 'shootingpointselector', 'volume', 'ensemble']:

            self.clone_storage(storage_name, storage2)

        storage2.close()


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
            storage_name = storage_to_copy.db

        for variable in self.variables.keys():
            if variable.startswith(storage_name + '_'):
                if variable not in new_storage.variables:
                    # orderparameters have additional variables in the storage that need to be copied
                    new_storage.createVariable(variable, str(self.variables[variable].dtype), self.variables[variable].dimensions)
                    for attr in self.variables[variable].ncattrs():
                        setattr(new_storage.variables[variable], attr, getattr(self.variables[variable], attr))

                    new_storage.variables[variable][:] = self.variables[variable][:]
                else:
                    for idx in range(0, len(self.variables[variable])):
                        new_storage.variables[variable][idx] = self.variables[variable][idx]



class StorableObjectJSON(paths.todict.ObjectJSON):
    def __init__(self, storage, unit_system = None, class_list = None):
        super(StorableObjectJSON, self).__init__(unit_system, class_list)
        self.excluded_keys = ['name', 'idx', 'json', 'identifier']
        self.storage = storage

    def simplify(self,obj, base_type = ''):
        if type(obj).__module__ != '__builtin__':
            if hasattr(obj, 'idx') and (not hasattr(obj, 'nestable') or (obj.base_cls_name != base_type)):
                # this also returns the base class name used for storage
                # store objects only if they are not creatable. If so they will only be created in their
                # top instance and we use the simplify from the super class ObjectJSON
                self.storage.save(obj)
                return { '_idx' : obj.idx[self.storage], '_cls' : obj.__class__.__name__ }

        return super(StorableObjectJSON, self).simplify(obj, base_type)

    def build(self,obj):
        if type(obj) is dict:
            if '_cls' in obj and '_idx' in obj:
                if obj['_cls'] in paths.todict.class_list:
                    base_cls = paths.todict.class_list[obj['_cls']]
                    result = self.storage.load(base_cls, obj['_idx'])
                else:
                    result = self.storage.load(obj['_cls'], obj['_idx'])

                return result

        return super(StorableObjectJSON, self).build(obj)