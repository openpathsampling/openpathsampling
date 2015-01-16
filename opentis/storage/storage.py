'''
Created on 06.07.2014

@author: JDC Chodera
@author: JH Prinz
'''

import netCDF4 as netcdf
import os.path

import numpy
import simtk.unit as u

import opentis as ops
import opentis.todict

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


    def _register_storages(self, store = None):
        if store is None:
            store = self

        self.trajectory = ops.storage.TrajectoryStorage(store).register()
        self.snapshot = ops.storage.SnapshotStorage(store).register()
        self.configuration = ops.storage.ConfigurationStorage(store).register()
        self.momentum = ops.storage.MomentumStorage(store).register()
        self.sample = ops.storage.SampleStorage(store).register()
        self.sampleset = ops.storage.SampleSetStorage(store).register()

        self.collectivevariable = ops.storage.ObjectDictStorage(store, ops.OrderParameter, ops.Configuration).register()

        self.pathmover = ops.storage.ObjectStorage(store, ops.PathMover, named=True, json=True, identifier='json').register()
        self.movedetails = ops.storage.ObjectStorage(store, ops.MoveDetails, named=False, json=True, identifier='json').register()
        self.shootingpoint = ops.storage.ObjectStorage(store, ops.ShootingPoint, named=False, json=True).register()
        self.shootingpointselector = ops.storage.ObjectStorage(store, ops.ShootingPointSelector, named=False, json=True, identifier='json').register()
        self.globalstate = ops.storage.ObjectStorage(store, ops.GlobalState, named=True, json=True, identifier='json').register()
        self.engine = ops.storage.ObjectStorage(store, ops.DynamicsEngine, named=True, json=True, identifier='json').register()
        self.volume = ops.storage.ObjectStorage(store, ops.Volume, named=True, json=True, identifier='json').register()
        self.ensemble = ops.storage.ObjectStorage(store, ops.Ensemble, named=True, json=True, identifier='json').register()

    def _setup_class(self):
        self._storages = {}
        self._storages_base_cls = {}
        self.links = []
        self.simplifier = ops.ObjectJSON()
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
                 template=None, n_atoms=None, units=None):
        '''
        Create a storage for complex objects in a netCDF file

        Parameters
        ----------
        filename : string
            filename of the netcdf file to be used or created
        mode : string, default: None
            the mode of file creation, one of 'w' (write), 'a' (append) or
            None, which will append any existing files.
        template : opentis.Snapshot
            a Snapshot instance that contains a reference to a Topology, the
            number of atoms and used units
        n_atoms : int or None
            If not None overrides the number of atoms in the storage
        units : dict of {str : simtk.unit.Unit } or None
            representing a dict of string representing a dimension ('length', 'velocity', 'energy') pointing to
            the simtk.unit.Unit to be used. If not None overrides the standard units used
        '''

        if mode == None:
            if os.path.isfile(filename):
                mode = 'a'
            else:
                mode = 'w'

        self.filename = filename

        super(Storage, self).__init__(filename, mode)

        self._setup_class()

        if units is not None:
            self.dimension_units.update(units)

        self._register_storages()

        if mode == 'w':
            self._initialize_netCDF()

            if template.topology is not None:
                self.topology = template.topology

            if n_atoms is not None:
                self.atoms = n_atoms
            elif self.topology is not None:
                self.atoms = self.topology.n_atoms
            elif template.coordinates is not None:
                self.atoms = template.coordinates.shape[0]
            else:
                raise RuntimeError("Storage given neither n_atoms nor topology")

            # update the units for dimensions from the template
            self.dimension_units.update(ops.tools.units_from_snapshot(template))
            self._init_storages(units=self.dimension_units)

            # create a json from the mdtraj.Topology() and store it
            self.write_str('topology', self.simplifier.to_json(self.simplifier.topology_to_dict(self.topology)))

            # Save the initial configuration
            self.snapshot.save(template)

            self.createVariable('template_idx', 'i4', 'scalar')
            self.variables['template_idx'][:] = template.idx[self]

            self.sync()

        elif mode == 'a' or mode == 'r+' or mode == 'r':
            self._restore_storages()

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

            self.topology = self.simplifier.topology_from_dict(self.simplifier.from_json(self.variables['topology'][0]))
            self.atoms = self.topology.n_atoms

    def __str__(self):
        return "OpenPathSampling netCDF Storage @ '" + self.filename + "'"

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
        return self.__dict__[item]

    def __setattr__(self, key, value):
        self.__dict__[key] = value

    def _init_storages(self, units=None):
        '''
        Run the initialization on all added classes, when the storage is created only!

        Notes
        -----
        Only runs when the storage is created.
        '''

        for storage in self.links:
            # create a member variable which is the associated Class itself
            storage.dimension_units.update(units)
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
            
        if 'spatial' not in self.dimensions:
            self.createDimension('spatial', 3) # number of spatial dimensions
        
        # Set global attributes.
        setattr(self, 'title', 'Open-Transition-Interface-Sampling')
        setattr(self, 'application', 'Host-Guest-System')
        setattr(self, 'program', 'run.py')
        setattr(self, 'programVersion', __version__)
        setattr(self, 'Conventions', 'Multi-State Transition Interface TPS')
        setattr(self, 'ConventionVersion', '0.1')

        # Create a string to hold the topology
        self.init_str('topology')
        self.write_str('topology', '')

        # Force sync to disk to avoid data loss.
        self.sync()

    def write_str(self, name, string):
        '''
        Write a string into a netCDF Variable

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
            the class name of the BaseClass of the stored object, which is needed when loading the object
            to identify the correct storage
        """

        if type(obj) is list:
            return [ self.save(part, *args, **kwargs) for part in obj]
        elif type(obj) is tuple:
            return tuple([self.save(part, *args, **kwargs) for part in obj])
        elif type(obj) in self._storages:
            store = self._storages[type(obj)]
            store.save(obj, *args, **kwargs)
            if type(obj) not in self._storages_base_cls:
                self._storages_base_cls[type(obj)] = obj.__class__.__name__
            return self._storages_base_cls[type(obj)]
        else:
            # Make sure subclassed objects will also be stored
            # Might come in handy someday
            for ty in self._storages:
                if type(ty) is not str and issubclass(type(obj), ty):
#                    print 'sub', ty, type(obj)
                    store = self._storages[ty]
                    # store the subclass in the _storages member for
                    # faster access next time
                    self._storages[type(obj)] = store
                    if type(ty) not in self._storages_base_cls:
                        self._storages_base_cls[ty] = ty.__name__
                    self._storages_base_cls[type(obj)] = self._storages_base_cls[ty]
                    store.save(obj, *args, **kwargs)

                    # this return the base_cls.__name__ to make sure on loading the right function is called
                    return self._storages_base_cls[ty]

        # Could not save this object. Might raise an exception, but return an empty string as type
        return ''

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
        If you want to load a subclassed Ensemble you need to load using `Ensemble` or `"Ensemble"`
        and not use the subclass
        """

        if obj_type in self._storages:
            store = self._storages[obj_type]
            return store.load(*args, **kwargs)
        elif type(obj_type) is type:
            # If we give a class we might be lucky and find the base class and
            # register is. If we try to load a type given by a string which is
            # not the base_class name originally registered there is no way to
            # tell.
            for ty in self._storages:
                if type(ty) is not str and issubclass(obj_type, ty):
                    store = self._storages[ty]
                    # store the subclass in the _storages member for
                    # faster access next time
                    self._storages[obj_type] = store
                    self._storages[obj_type.__name__] = store

                    if type(ty) not in self._storages_base_cls:
                        self._storages_base_cls[ty] = ty.__name__
                        self._storages_base_cls[ty.__name__] = ty.__name__

                    self._storages_base_cls[obj_type] = self._storages_base_cls[ty]
                    self._storages_base_cls[obj_type.__name__] = ty.__name__

                    store = self._storages[obj_type]
                    return store.load(*args, **kwargs)

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


    def clone_storage(self, storage_to_copy, new_storage):
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



class StorableObjectJSON(ops.todict.ObjectJSON):
    def __init__(self, storage, unit_system = None, class_list = None):
        super(StorableObjectJSON, self).__init__(unit_system, class_list)
        self.excluded_keys = ['name', 'idx', 'json', 'identifier']
        self.storage = storage

    def simplify(self,obj, base_type = ''):
        if type(obj).__module__ != '__builtin__':
#            print obj.__dict__, hasattr(obj, 'creatable')
            if hasattr(obj, 'idx') and (not hasattr(obj, 'nestable') or (obj.base_cls_name != base_type)):
                # this also returns the base class name used for storage
                # store objects only if they are not creatable. If so they will only be created in their
                # top instance and we use the simplify from the super class ObjectJSON
                base_cls = self.storage.save(obj)
                return { '_idx' : obj.idx[self.storage], '_base' : base_cls, '_cls' : obj.__class__.__name__ }

        return super(StorableObjectJSON, self).simplify(obj, base_type)

    def build(self,obj):
        if type(obj) is dict:
            if '_base' in obj and '_idx' in obj:
                result = self.storage.load(obj['_base'], obj['_idx'])
                # restore also the actual class name

                result.cls = obj['_cls']
                return result

        return super(StorableObjectJSON, self).build(obj)