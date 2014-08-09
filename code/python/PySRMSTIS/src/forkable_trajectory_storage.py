'''
Created on 06.07.2014

@author: JDC Chodera
@author: JH Prinz
'''

import numpy

import netCDF4 as netcdf # for netcdf interface provided by netCDF4 in enthought python

import simtk.unit as units

from simtk.unit import nanosecond, picosecond, nanometers, nanometer, picoseconds, femtoseconds, femtosecond, kilojoules_per_mole, Quantity

from trajectory import Trajectory
from snapshot import Snapshot

import pickle
import mdtraj as md
import inspect

import os.path

# TODO: Remove all stuff that is content related and allow to register a class with the storage

#=============================================================================================
# SOURCE CONTROL
#=============================================================================================

__version__ = "$Id: NoName.py 1 2014-07-06 07:47:29Z jprinz $"



#=============================================================================================
# Multi-State Transition Interface Sampling
#=============================================================================================

class ForkableTrajectoryStorage(object):
    '''
    A netCDF4 wrapper to store trajectories based on snapshots of an OpenMM simulation. This allows effective storage of shooting trajectories
    '''


    def __init__(self, topology, filename = 'trajectory.nc', mode = 'auto'):
        '''
        Create a storage for trajectories that are stored by snapshots and list of their indices
        
        ARGUMENTS
        
        topology (openmm.app.Topology) - the topology of the system to be stored. Needed for 
        filename (string) - filename of the netcdf file
        
        OPTIONAL ARGUMENTS
        
        mode (string, default: 'auto') - the mode of file creation, one of 'create', 'restore' or 'auto'
        '''
        
        self.fn_storage = filename

        # if not specified then create if not existing
        if mode == 'auto':                    
            if os.path.isfile(filename):
                mode = 'restore'
            else:
                mode = 'create'
        
        if mode == 'create':
            self.topology = md.Topology.from_openmm(topology)
            self._initialize_netcdf()
        elif mode == 'restore':
            self._restore_netcdf()
            
        self.links = []
        
    def add(self, class_obj):
        
        if hasattr(class_obj, '_initialize_netCDF') and hasattr(class_obj, '_restore_netCDF'):
            #class is compatible and has necessary classes
            self.links.append(class_obj)
            
        pass
    
    def init_classes(self):
        for cls in self.links:
            cls._initialize_netcdf(self)
            
            # add all static methods fromt cls to self.'cls'.method_name
            
#            methods = inspect.getmembers(cls, predicate)

    def _initialize_netcdf(self):
        """
        Initialize the netCDF file for storage.
        
        """    

        # Open NetCDF file for writing
        ncfile = netcdf.Dataset(self.fn_storage, 'w') # for netCDF4
        
        # Store netcdf file handle.
        self.ncfile = ncfile

        if 'scalar' not in self.ncfile.dimensions:
            self.ncfile.createDimension('scalar', 1) # scalar dimension
        
        # Set global attributes.
        setattr(ncfile, 'title', 'Multi-State-Transition-Interface-Sampling')
        setattr(ncfile, 'application', 'Host-Guest-System')
        setattr(ncfile, 'program', 'run.py')
        setattr(ncfile, 'programVersion', __version__)
        setattr(ncfile, 'Conventions', 'Multi-State Transition Interface TPS')
        setattr(ncfile, 'ConventionVersion', '0.1')
        
        # initialize arrays used for snapshots
        Snapshot._init_netcdf(self)

        # initialize arrays used for trajectories
        Trajectory._init_netcdf(self)
                        
        # Force sync to disk to avoid data loss.
        ncfile.sync()

        return
    
    def _restore_netcdf(self):
        '''
        Restore the storage from the netCDF file
        '''
        # Open NetCDF file for appending
        ncfile = netcdf.Dataset(self.fn_storage, 'a')
        
        # Store netcdf file handle.
        self.ncfile = ncfile
        
        # initialize arrays used for snapshots
        Snapshot._restore_netcdf(self)

        # initialize arrays used for trajectories
        Trajectory._restore_netcdf(self)
        
        return
    
    def snapshot(self, idx):
        return Snapshot.load(idx)
    
    def snapshot_coordinates_as_array(self, frame_indices=None, atom_indices=None):
        return Snapshot.coordinates_as_numpy(self, frame_indices, atom_indices)
        
    def trajectory_coordinates_as_array(self, idx, atom_indices=None):
        frame_indices = Trajectory.load_indices(idx)
        
        return self.snapshot_coordinates_as_array(frame_indices, atom_indices)        
    
    def trajectory(self, idx):        
        return Trajectory.load(idx)
    
    def number_of_trajectories(self):
        return Trajectory.load_number()

    def number_of_snapshots(self):
        return Snapshot.load_number()
    
    def last_trajectory(self):
        return self.trajectory(self.number_of_trajectories())

    def all_trajectory_indices(self):
        '''
        Return a list of list of frame indices
        '''
        return Trajectory.load_all_indices()
    
    def all_snapshot_coordinates_as_mdtraj(self, atom_indices = None):
        """
        Return all snapshots as a mdtraj.Trajectory object using only the specified atoms
        
        OPTIONAL ARGUMENTS
        
        atom_indices (list of int, Default: None) - list of atom indices to be used for the trajectory
         
        """
                    
        output = self.snapshot_coordinates_as_array(atom_indices = atom_indices)
                
        topology = self.topology
                
        if atom_indices is not None:
            topology = topology.subset(atom_indices)
                                                 
        return md.Trajectory(output, topology)
        
    def _store_metadata(self, groupname, metadata):
        """
        Store metadata in NetCDF file.
        
        """

        # Create group.
        ncgrp = self.ncfile.createGroup(groupname)

        # Store metadata.
        for (key, value) in metadata.iteritems():
            # TODO: Handle more sophisticated types.
            ncvar = ncgrp.createVariable(key, type(value))
            ncvar.assignValue(value)

        return
        
    def _store_options(self, obj, group_name = 'options'):
        """
        Store run parameters in NetCDF file.

        """

        self.verbose_root = False

        # Create a group to store state information.
        ncgrp_options = self.ncfile.createGroup(group_name)
        
        # Store run parameters.
        for option_name in obj.options_to_store:
            # Get option value.
            option_value = getattr(obj, option_name)            
            self._store_single_option(ncgrp_options, option_name, option_value)
            
            if self.verbose_root: print "Storing option: %s -> %s (type: %s)" % (option_name, option_value, type(option_value))            

        return

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
            ncvar = ncgrp.createVariable(obj_name, type(obj_value))
            ncvar.assignValue(obj_value)
            
        if option_unit: 
            setattr(ncvar, 'units', str(option_unit))
 
        setattr(ncvar, 'type', option_type.__name__) 
        setattr(ncvar, 'module', option_type.__module__)         

        return


    def _restore_options(self, obj, group_name = 'options'):
        """
        Restore run parameters from NetCDF file.
        """
        
        self.verbose_root = False

        if self.verbose_root: print "Attempting to restore options from NetCDF file..."

        # Make sure this NetCDF file contains option information
        if not group_name in self.ncfile.groups:
            # Not found, signal failure.
            return False

        # Find the group.
        ncgrp_options = self.ncfile.groups[group_name]
        
        print ncgrp_options.variables.keys()

        # Load run parameters.
        for option_name in ncgrp_options.variables.keys():
            # Get NetCDF variable.
            option_value = self._restore_single_option(ncgrp_options, option_name)
            
            print option_name
            
            # Store option.
            if self.verbose_root: print "Restoring option: %s -> %s (type: %s)" % (option_name, str(option_value), type(option_value))
            setattr(obj, option_name, option_value)
            
        # Signal success.
        return True
    
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