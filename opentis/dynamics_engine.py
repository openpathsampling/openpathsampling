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
from opentis.snapshot import Configuration, Momentum
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
    Class to wrap a simulation tool to store the context and rerun, needed
    parameters, storage, etc. 
    
    Notes
    -----
    Should be considered an abstract class: only its subclasses can be
    instantiated.
    '''


    def __init__(self, filename=None, opts=None, mode='auto'):
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

        # set up the opts
        self.opts = {}
        self.add_stored_parameters(opts)

        if not hasattr(self, 'topology'):
            self.topology = None

        # storage
        self.fn_storage = filename 

        # tell everybody who their engine is
        Snapshot.engine = self
        Configuration.engine = self
        Trajectory.engine = self

        if mode == 'create':
            # set up the max_length_stopper (if n_frames_max is given)
            # TODO: switch this not needing slice; use can_append
            if hasattr(self, 'n_frames_max'):
                self.max_length_stopper = LengthEnsemble(slice(0,self.n_frames_max-1))

            # storage
            self.storage = Storage(
                topology_file=self.topology,
                filename=self.fn_storage,
                mode='w'
            )
            self.storage.engine = self
            self.storage._store_options(self)
            Trajectory.storage = self.storage
 

    def add_stored_parameters(self, param_dict):
        '''Adds parameters in param_dict to the attribute dictionary of the
        DynamicsEngine object, and saves the relevant keys as options to
        store.

        Parameters
        ----------
        param_dict : dict
            dictionary of attributes to be added (and stored); attribute
            names are keys, with appropriate values
        '''
        if param_dict is not None:
            for key in param_dict.keys():
                self.opts[key] = param_dict[key]
                setattr(self, key, param_dict[key])
            self.options_to_store = self.opts.keys()
        return
 
    def start(self, snapshot=None):
        if snapshot is not None:
            self.current_snapshot = snapshot

    def stop(self, trajectory):
        """Nothing special needs to be done for direct-control simulations
        when you hit a stop condition."""
        pass

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
            self.start()

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

