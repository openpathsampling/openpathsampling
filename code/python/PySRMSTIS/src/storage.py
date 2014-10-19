'''
Created on 06.07.2014

@author: JDC Chodera
@author: JH Prinz
'''

import netCDF4 as netcdf # for netcdf interface provided by netCDF4 in enthought python
import pickle
import os.path

import numpy
import simtk.unit as units

from trajectory_store import TrajectoryStorage
from snapshot_store import SnapshotStorage, ConfigurationStorage, MomentumStorage

import mdtraj as md


#=============================================================================================
# SOURCE CONTROL
#=============================================================================================

__version__ = "$Id: NoName.py 1 2014-07-06 07:47:29Z jprinz $"

#=============================================================================================
# NetCDF Storage for multiple forked trajectories
#=============================================================================================

class Storage(netcdf.Dataset):
    '''
    A netCDF4 wrapper to store trajectories based on snapshots of an OpenMM simulation. This allows effective storage of shooting trajectories
    '''

    def __init__(self, filename = 'trajectory.nc', mode = None, topology = None):
        '''
        Create a storage for complex objects in a netCDF file
        
        Parameters
        ----------        
        topology : openmm.app.Topology
            the topology of the system to be stored. Needed for 
        filename : string
            filename of the netcdf file
        mode : string, default: None
            the mode of file creation, one of 'w' (write), 'a' (append) or None, which will append any existing files.
        '''

        if mode == None:
            if os.path.isfile(filename):
                mode = 'a'
            else:
                mode = 'w'

        self.filename = filename

        super(Storage, self).__init__(filename, mode)
        self._register_classes()

        if mode == 'w':
            self._init()

            if type(topology) is str:
                self.topology = md.load(topology)
            else:
                self.topology = topology

            self.atoms = self.topology.n_atoms
            self.topology.save_pdb('tempXXX.pdb')

            with open ('tempXXX.pdb', "r") as myfile:
                pdb_string=myfile.read()

            self.variables['pdb'][0] = pdb_string
            self._init_classes()
            self.sync()

        elif mode == 'a':
            self.pdb = self.variables['pdb'][0]

            with open ('tempXXX.pdb', "w") as myfile:
                myfile.write(self.pdb)

            self.topology = md.load('tempXXX.pdb')
            self._restore_classes()


    def __getattr__(self, item):
        return self.__dict__[item]

    def __setattr__(self, key, value):
        self.__dict__[key] = value


    def _add_class(self, class_obj):
        '''
        Add a class to the storage
        
        Parameters
        ----------
        class_obj : Class
            the class to be added
        '''

        if hasattr(class_obj, '_init'):
            self.links.append(class_obj)

    def _register_classes(self):
        '''
        Register all class types to be storable in the netCDF file.

        Notes
        -----
        If you want more, add it here!
        '''

        self.links = []
        self._add_class(MomentumStorage)
        self._add_class(ConfigurationStorage)
        self._add_class(SnapshotStorage)
        self._add_class(TrajectoryStorage)

        for cls in self.links:
            # create a member variable which is the associated Class itself
            store_name = cls.__name__[:-7].lower()

            setattr(self, store_name, cls(self))
            getattr(self, store_name).storage = self

    def _init_classes(self):
        '''
        Run the initialization on all added classes, when the storage is created only!

        Notes
        -----
        Only runs when the storage is created.
        '''

        for cls in self.links:
            # create a member variable which is the associated Class itself
            store_name = cls.__name__[:-7].lower()
            getattr(self, store_name)._init()

    def _restore_classes(self):
        '''
        Run restore on all added classes. Usually there is nothing to do.
        '''
        for cls in self.links:
            store_name = cls.__name__[:-7].lower()
#            getattr(self, store_name)._restore()


    def _init(self):
        """
        Initialize the netCDF file for storage itself.
        """

        # TODO: Allow to set the project parameters somehow!

        # _add_class shared dimension for everyone. scalar and spatial
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

        ncvar = self.createVariable('pdb', 'str', 'scalar')
        packed_data = numpy.empty(1, 'O')
        packed_data[0] = ""
        ncvar[:] = packed_data
                        
        # Force sync to disk to avoid data loss.
        self.sync()

    def _store_metadata(self, groupname, metadata):
        """
        Store metadata in NetCDF file.
        
        """

        # Create group.
        ncgrp = self.createGroup(groupname)

        # Store metadata.
        for (key, value) in metadata.iteritems():
            # TODO: Handle more sophisticated types.
            ncvar = ncgrp.createVariable(key, type(value))
            ncvar.assignValue(value)


    def _store_options(self, obj, group_name = 'options'):
        """
        Store run parameters in NetCDF file.

        """

        self.verbose_root = False

        # Create a group to store state information.
        ncgrp_options = self.createGroup(group_name)
        
        # Store run parameters.
        for option_name in obj.options_to_store:
            # Get option value.
            option_value = getattr(obj, option_name)            
            self._store_single_option(ncgrp_options, option_name, option_value)
            
            if self.verbose_root: print "Storing option: %s -> %s (type: %s)" % (option_name, option_value, type(option_value))            

    def _restore_options(self, obj, group_name = 'options'):
        """
        Restore run parameters from NetCDF file.
        """
        
        self.verbose_root = False

        if self.verbose_root: print "Attempting to restore options from NetCDF file..."

        # Make sure this NetCDF file contains option information
        if not group_name in self.groups:
            # Not found, signal failure.
            return False

        # Find the group.
        ncgrp_options = self.groups[group_name]

        # Load run parameters.
        for option_name in ncgrp_options.variables.keys():
            # Get NetCDF variable.
            option_value = self._restore_single_option(ncgrp_options, option_name)
            
            # Store option.
            if self.verbose_root: print "Restoring option: %s -> %s (type: %s)" % (option_name, str(option_value), type(option_value))
            setattr(obj, option_name, option_value)
            
        # Signal success.
        return True

    def _store_single_option(self, ncgrp, obj_name, obj_value):
        """
        Store run parameters in NetCDF file.

        """
                
        # If Quantity, strip off units first.
        option_unit = None
        if type(obj_value) == units.Quantity:
            option_unit = obj_value.unit
            obj_value = obj_value / option_unit
        # Store the Python type.
        option_type = type(obj_value)
        # Handle booleans
        if type(obj_value) == bool:
            obj_value = int(obj_value)
        # Store the variable.


        if option_type.__module__ != '__builtin__':
            # non of the standard cases -> use pickle to serialize
            obj_value = pickle.dumps(obj_value)
            
        if option_type is list:
            obj_value = pickle.dumps(obj_value)

        if option_type is dict:
            obj_value = pickle.dumps(obj_value)

        if option_type is tuple:
            obj_value = pickle.dumps(obj_value)
            
            
        if type(obj_value) == str:
            ncvar = ncgrp.createVariable(obj_name, type(obj_value), 'scalar')
            packed_data = numpy.empty(1, 'O')
            packed_data[0] = obj_value
            ncvar[:] = packed_data
        else:
            ncvar = ncgrp.createVariable(obj_name, type(obj_value), 'scalar')
            packed_data = numpy.empty(1, 'O')
            packed_data[0] = obj_value
            ncvar[:] = packed_data

        if option_unit: 
            setattr(ncvar, 'units', str(option_unit))
 
        setattr(ncvar, 'type', option_type.__name__) 
        setattr(ncvar, 'module', option_type.__module__)         

        return

    
    def _restore_single_option(self, ncgrp, obj_name):
        """
        Restore run parameters from NetCDF file.
        """
        
        # Get NetCDF variable.
        option_ncvar = ncgrp.variables[obj_name]
        # Get option value.
        obj_value = option_ncvar[0]
        # Cast to Python type.
        type_name = getattr(option_ncvar, 'type') 
        module_name = getattr(option_ncvar, 'module')         
        if module_name != '__builtin__':
            obj_value = str(obj_value)
            obj_value = pickle.loads(obj_value)
        elif type_name == 'tuple' or type_name == 'dict' or type_name == 'list':
            obj_value = str(obj_value)
            obj_value = pickle.loads(obj_value)
        else:
            obj_value = eval(type_name + '(' + repr(obj_value) + ')')
            # If Quantity, assign units.
            if hasattr(option_ncvar, 'units'):
                option_unit_name = getattr(option_ncvar, 'units')
                if option_unit_name[0] == '/': 
                    option_unit_name = '1' + option_unit_name
                option_unit = eval(option_unit_name, vars(units))
                obj_value = units.Quantity(obj_value, option_unit)

        # Signal success.
        return obj_value