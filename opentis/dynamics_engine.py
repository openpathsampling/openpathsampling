'''
Created on 01.07.2014

@author JDC Chodera
@author: JH Prinz
'''

import time

import numpy as np
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import mdtraj as md

from opentis.storage import Storage
from opentis.trajectory import Trajectory
from opentis.snapshot import Snapshot
from opentis.integrators import VVVRIntegrator
from opentis.ensemble import LengthEnsemble
#from opentis.openmm_simulation import OpenMMSimulation


#=============================================================================
# SOURCE CONTROL
#=============================================================================

__version__ = "$Id: NoName.py 1 2014-07-06 07:47:29Z jprinz $"

#=============================================================================
# Multi-State Transition Interface Sampling
#=============================================================================

class DynamicsEngine(object):
    '''
    Class to wrap a simulation tool to store the context and rerun, needed parameters, storage, etc. 
    
    Notes
    -----
    - Should only be contructed through factory functions
    - This is the main object knowing about all things of the simulation. It
      might be possible to replace this with other simulation tools than
      OpenMM
    
    '''


    def __init__(self):
        '''
        Create an empty DynamicsEngine object
        
        Notes
        -----
        The purpose of an engine is to create trajectories and keep track
        of the results. The main method is 'generate' to create a
        trajectory, which is a list of snapshots and then can store the in
        the associated storage. In the initialization this storage is
        created as well as the related Trajectory and Snapshot classes are
        initialized.
        '''
        self.op = None
        self.storage = None
        self.initialized = False
        self.running = dict()

    def generate(self, snapshot, running = None):
        r"""
        Generate a velocity Verlet trajectory consisting of ntau segments of
        tau_steps in between storage of Snapshots and randomization of
        velocities.

        Parameters
        ----------
        snapshot : Snapshot 
            initial coordinates; velocities will be assigned from
            Maxwell-Boltzmann distribution            
        running : list of function(Snapshot)
            callable function of a 'Snapshot' that returns True or False.
            If one of these returns False the simulation is stopped.

        Returns
        -------    
        trajectory : Trajectory
            generated trajectory of initial conditions, including initial
            coordinate set

        Notes
        -----
        Might add a return variable of the reason why the trajectory was
        aborted. Otherwise check the length and compare to max_frames
        """

        # Are we ready to rumble ?
        if self.initialized:
            
            self.current_snapshot = snapshot

            # Store initial state for each trajectory segment in trajectory.
            trajectory = Trajectory()
            trajectory.append(snapshot)
            
            frame = 0
            
            stop = False

            while stop == False:
                                
                # Do integrator x steps
                snapshot = self.generate_next_frame()
                frame += 1
                
                # Store snapshot and add it to the trajectory. Stores also
                # final frame the last time
                if self.storage is not None:
                    self.storage.snapshot.save(snapshot)
                trajectory.append(snapshot)
                
                # Check if we should stop. If not, continue simulation
                if running is not None:
                    for runner in running:
                        #print str(runner), runner(trajectory)
                        keep_running = runner(trajectory)
                        self.running[runner] = keep_running
                        stop = stop or not keep_running

                # TODO: switch to self.max_length_stopper.can_append()
                stop = stop or not self.max_length_stopper(trajectory)
                
            self.stop(trajectory)
            return trajectory
        else:
            raise RuntimeWarning("Can't generate from an uninitialized system!")
            return None

