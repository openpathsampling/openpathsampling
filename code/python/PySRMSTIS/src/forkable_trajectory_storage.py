'''
Created on 06.07.2014

@author: JDC Chodera
@author: JH Prinz
'''

import numpy

import netCDF4 as netcdf # for netcdf interface provided by netCDF4 in enthought python

from simtk.unit import nanosecond, picosecond, nanometers, nanometer, picoseconds, femtoseconds, femtosecond, kilojoules_per_mole, Quantity

from trajectory import Trajectory
from snapshot import Snapshot

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


    def __init__(self, params):
        '''
        Constructor
        '''
        
        self.store_filename = "database.nc"
        self.context = None
        
        self._initialize_netcdf()
        
    def _initialize_netcdf(self):
        """
        Initialize NetCDF file for storage.
        
        """    

        # Open NetCDF file for writing
        ncfile = netcdf.Dataset(self.store_filename, 'w') # for netCDF4
        
        # Store netcdf file handle.
        self.ncfile = ncfile

        # Set global attributes.
        setattr(ncfile, 'tile', self.title)
        setattr(ncfile, 'application', 'Host-Guest-System')
        setattr(ncfile, 'program', 'run.py')
        setattr(ncfile, 'programVersion', __version__)
        setattr(ncfile, 'Conventions', 'Multi-State Transition Interface TPS')
        setattr(ncfile, 'ConventionVersion', '0.1')

        # TODO: save information about the system to be simulated for restart... ???
        setattr(ncfile, 'temperature', 300.0); # in Kelvin

        
        # initialize arrays used for snapshots
        Snapshot._init_netcdf(ncfile, self.context)

        # initialize arrays used for trajectories
        Trajectory._init_netcdf(ncfile)
                        
        # Force sync to disk to avoid data loss.
        ncfile.sync()

        return
    
    def snapshot(self, idx):
                
        x = self.ncfile.variables['trajectory_coordinates'][idx,:,:].astype(numpy.float32).copy()
        coordinates = Quantity(x, nanometers)                
        v = self.ncfile.variables['trajectory_velocities'][idx,:,:].astype(numpy.float32).copy()
        velocities = Quantity(v, nanometers / picoseconds)                
        V = self.ncfile.variables['trajectory_potential'][idx]
        potential_energy = Quantity(V, kilojoules_per_mole)
        T = self.ncfile.variables['trajectory_kinetic'][idx]
        kinetic_energy = Quantity(T, kilojoules_per_mole)
        snapshot = Snapshot(coordinates=coordinates, velocities=velocities, kinetic_energy=kinetic_energy, potential_energy=potential_energy)
        snapshot.idx = idx
        return snapshot
    
    def trajectory(self, idx):        
        
        nframes = self.ncfile.variables['trajectory_length'][idx]

        trajectory = Trajectory()
        for frame_index in range(nframes):                
            snapshot = self.snapshot( self.ncfile.variables['trajectory_idx'][frame_index])
            trajectory.append(snapshot)
            
        return trajectory
    
    
    def _resume_from_netcdf(self):
        """
        Resume execution by reading current positions and energies from a NetCDF file.
        
        """
        
        # NOT IMPLEMENTED YET!
        return

        # Open NetCDF file for reading
        ncfile = netcdf.Dataset(self.store_filename, 'r') # for netCDF4
        
        # TODO: Perform sanity check on file before resuming

        # Get current dimensions.
        self.iteration = ncfile.variables['activities'].shape[0] - 1
        self.nstates = ncfile.variables['activities'].shape[1]
        self.nframes = ncfile.variables['trajectory_coordinates'].shape[1]
        self.natoms = ncfile.variables['trajectory_coordinates'].shape[2]

        print "iteration = %d, nstates = %d, natoms = %d" % (self.iteration, self.nstates, self.natoms)

        # Close NetCDF file.
        ncfile.close()        
        
        # Reopen NetCDF file for appending, and maintain handle.
        self.ncfile = netcdf.Dataset(self.store_filename, 'a') # for netCDF4
        
        
        return