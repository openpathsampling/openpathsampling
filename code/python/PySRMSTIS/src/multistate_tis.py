'''
Created on 01.07.2014

@author JDC Chodera
@author: JH Prinz
'''


import numpy

import netCDF4 as netcdf # for netcdf interface provided by netCDF4 in enthought python


from simtk.unit import nanosecond, picosecond, nanometers, nanometer, picoseconds, femtoseconds, femtosecond, kilojoules_per_mole, Quantity

from trajectory import Trajectory
from snapshot import Snapshot
from forkable_trajectory_storage import ForkableTrajectoryStorage


#=============================================================================================
# SOURCE CONTROL
#=============================================================================================

__version__ = "$Id: NoName.py 1 2014-07-06 07:47:29Z jprinz $"



#=============================================================================================
# Multi-State Transition Interface Sampling
#=============================================================================================


class MultiStateTIS(object):
    '''
    Class to run multistate transition interface sampling simulations
    
    LITERATURE
    
    NOTES
    
    Needs to know a way to define coresets and surrounding regions (regions between interfaces). This definition should be kept
    seperate from the algorithm to be flexible. We want to use MSM. So the general workflow is
    
    1. generate TIS trajectories that go between states and try to leave states and sample all known interfaces equally likely.
    2. analyze all generated trajectories to decide if new cores have appeared and if we need to adjust the interface regions
    3. go to 1 until we are satisfied with the generated model 
    
    If want to utilize MSMs in the process for step 2. Once we have a MSM we can use it to compute lots of key properties.
    
    What we need.
    - a way to generate new trajectories (from TPS)
    - a way to cluster all data into a MSM and get back the state definition.
    
    Class Objects needed.
    - trajectory: list of snapshots        [ok]
    - snapshot                             [ok]
    - storage using netCDF                 [ok]
    - TPS with stopping                    
    - Clustering
    - MSM
    - MSM Analysis
    - region predicter
    - 
    '''


    def __init__(self, params):
        '''
        Constructor
        
        '''
        
        self.xframes = 1000;            # maximal length of trajectories in saved frames, needed for 
        self.natoms = 10;               # number of atoms in the system
        
        Snapshot.context = self.context
        
        # Create a trajectory storage
        # TODO: Add naming capabilities, other options and add the system to be simulated
        self.storage = ForkableTrajectoryStorage()
        