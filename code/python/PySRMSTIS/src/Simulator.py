'''
Created on 01.07.2014

@author JDC Chodera
@author: JH Prinz
'''


import numpy as np
import time

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *

import mdtraj as md


from simtk.unit import nanosecond, picosecond, nanometers, nanometer, picoseconds, femtoseconds, femtosecond, kilojoules_per_mole, Quantity

from trajectory import Trajectory
from snapshot import Snapshot
from storage import Storage
from integrators import VVVRIntegrator
from ensemble import LengthEnsemble

#=============================================================================================
# SOURCE CONTROL
#=============================================================================================

__version__ = "$Id: NoName.py 1 2014-07-06 07:47:29Z jprinz $"

#=============================================================================================
# Multi-State Transition Interface Sampling
#=============================================================================================

class Simulator(object):
    '''
    Class to wrap a simulation tool to store the context and rerun, needed parameters, storage, etc. 
    
    Notes
    -----
    
    - Should only be contructed through factory functions
    - This is the main object knowing about all things of the simulation. It might be possible to replace this with other simulation tools than OpenMM
    
    '''


    def __init__(self):
        '''
        Create an empty simulator object
        
        Notes
        -----
        The purpose of a simulator is to create simulations and keep track of the results. The main method is 'generate' to create a trajectory, which is
        a list of snapshots and then can store the in the associated storage. In the initialization this storage is created as well as the related
        Trajectory and Snapshot classes are initialized.
        
        '''
        self.op = None
        self.initialized = False
        
    def generate(self, snapshot, running = None):
        r"""
        Generate a velocity Verlet trajectory consisting of ntau segments of tau_steps in between storage of Snapshots and randomization of velocities.

        Parameters
        ----------
        snapshot : Snapshot 
            initial coordinates; velocities will be assigned from Maxwell-Boltzmann distribution            
        running : function(Snapshot) 
            callable function of a 'Snapshot' that returns True or False

        Returns
        -------    
        trajectory : Trajectory
            generated trajectory of initial conditions, including initial coordinate set

        Notes
        -----
        Might add a return variable of the reason why the trajectory was aborted. Otherwise check the length and compare to max_frames
        """

        # Are we ready to rumble ?
        if self.initialized:
            # Store initial state for each trajectory segment in trajectory.
            trajectory = Trajectory()
            
            # Set initial positions
            self.simulation.context.setPositions(snapshot.coordinates)

            trajectory.append(snapshot)
            
            # Assign velocities from Maxwell-Boltzmann distribution          
            self.simulation.context.setVelocities(snapshot.velocities)
#            self.simulation.context.setVelocitiesToTemperature(self.temperature)
        
            # Propagate dynamics by velocity Verlet.            
            frame = 0
            nsteps_per_iteration = self.nframes_per_iteration            
            
            stop = False

            while stop == False:
                                
                # Do integrator x steps
                self.simulation.step(nsteps_per_iteration)            
                frame += 1
                
                # Store snapshot and add it to the trajectory. Stores also final frame the last time
                snapshot = Snapshot(self.simulation.context)
                self.storage.snapshot.save(snapshot)
                trajectory.append(snapshot)
                
                # Check if reached a core set. If not, continue simulation
                if running is not None:
                    for runner in running:
#                        print str(runner), runner(trajectory)
                        stop = stop or not runner(trajectory)

                # We could also just count the number of frames. Might be faster but not as nice :)                
                stop = stop or not self.max_length_stopper(trajectory)
                
                if self.op:
                    print self.max_length_stopper, self.max_length_stopper(trajectory)
                    print frame
                    print len(trajectory)
                    print [ s.idx for s in trajectory]
                    
                    print 'OP :', self.op(snapshot)
                
                    
            return trajectory
        else:
            # TODO: Throw an error! Needs to be initialized
            return None

    @staticmethod
    def Alanine_system(mode ='auto'):
        '''
        Setup / Restore a simulator for the Alanine Dipeptide system
        
        OPTIONAL ARGUMENTS
        
        mode (string) - Set to 'restore' to restore from file, 'create' to create fresh simulation and 'auto' to restore only if no database is present
        '''
        
        self = Simulator()
        self.fn_storage = "data/trajectory.nc"
        
        # associate Snapshots and Trajecotories with this simulation
        Snapshot.simulator = self
        Trajectory.simulator = self
        
        if mode == 'auto':
            if os.path.isfile(self.fn_storage):
                mode = 'restore'
            else:
                mode = 'create'

        if mode == 'create':
            self._set_alanine_options()
            self._create_OpenMMSimulation()
            
#            self.system_serial = openmm.XmlSerializer.serialize(self.system)
            
            self._equilibrate_system()
            
            # Create a trajectory storage
            self.storage = Storage(
                                                 topology = self.fn_initial_pdb,
                                                 filename = self.fn_storage,
                                                 mode = 'w'
                                                 )
            self.storage.simulator = self
            # save options
            self.storage._store_options(self)
            
            # save initial equilibrated frame as snapshot ID #0. Might be useful later, who knows
            snapshot = Snapshot(self.simulation.context)

            self.storage.snapshot.save(snapshot, 0,0)
        
        if mode == 'restore':
            # Need the oposite order, first open database 
            self.storage = Storage(
                                                     topology = None,
                                                     filename = self.fn_storage,
                                                     mode = 'a'
                                                     )

            # and load options
            self.storage._restore_options(self)        
            
            # This is not stored yet so we need to set it again. Should be stored in the database along with the still missing topology
#            self.solute_indices = range(22)
            self.platform = 'CPU'

            # Finally create the OpenMM simulation system
            self._create_OpenMMSimulation()
                        
            # Still missing is the topology which is still not in the database and can only be accessed after the simulation has been created
            self.storage.topology = md.Topology.from_openmm(self.topology)
            
        # Finished
        self.initialized = True
        
        return self
            
    def _set_alanine_options(self):        
        self.options_to_store = [ 
                                 'temperature', 
                                 'collision_rate', 
                                 'timestep', 
                                 'forcefield_solute', 
                                 'forcefield_solvent', 
                                 'nframes_per_iteration',
                                 'n_frames_max',
                                 'start_time',
                                 'fn_initial_pdb',
                                 'integrator_class',
                                 'fn_storage',
                                 'platform',
                                 'topology',
                                 'solute_indices',
                                 'system_serial',
                                 'integrator_serial',
                                 'pdb_file'
                                 ]
        
        self.temperature = 300.0 * kelvin                       # temperature
        self.collision_rate = 1.0 / picoseconds                 # collision rate for VVVR integrator
        self.timestep = 2.0 * femtoseconds                      # timestep for Langevin integrator
        
        self.nframes_per_iteration = 10                         # number of steps of dynamics per iteration
        
        self.fn_initial_pdb = 'data/Alanine_solvated.pdb'
        self.forcefield_solute = 'amber99sbildn.xml'
        self.forcefield_solvent = 'tip3p.xml'        
        self.platform = 'CUDA'
        self.platform = 'CPU'
        self.solute_indices = range(22)                         # indices of Alanine without water

        self.n_frames_max = 5000;                               # maximal length of trajectories in saved frames
        
        self.start_time = time.time()                           # the time when we started
        
        self.pdb_file = PDBFile(self.fn_initial_pdb)
        self.topology = self.pdb_file.topology


    def _create_OpenMMSimulation(self):
        '''
        Create an OpenMM simulation from the parameters specified before.
        
        NOTES
        
        This should be generic and not depend on other parameters than then ones stored in the netcdf file. Still missing is the integrator type, the topology and the solute coordinates
        which at some point might be accessible from the topology
        '''

        platform = openmm.Platform.getPlatformByName(self.platform)    # platform to use
        
        #TODO: Serialize initial pdb or better the used topology so that we do not need to read an external file
        forcefield = ForceField(self.forcefield_solute, self.forcefield_solvent)
        
        system = forcefield.createSystem(self.topology, nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)
        
        
        integrator = VVVRIntegrator(self.temperature, self.collision_rate, self.timestep)
        self.integrator_class = type(integrator).__name__
        
        simulation = Simulation(self.topology, system, integrator, platform)
        
        self.system_serial = openmm.XmlSerializer.serialize(system)
        self.integrator_serial = openmm.XmlSerializer.serialize(integrator)
        
#        print self.system_serial
                
        # This could be moved to initialization since it will not change
        self.max_length_stopper = LengthEnsemble(slice(0,self.n_frames_max - 1))        
                                
        self.simulation = simulation
                
    def _equilibrate_system(self):
        '''
        Equilibrate the Alanine System. Just a helper to keep the code clean. One could also provide an equilibrated structure from the start
        '''

        #=============================================================================================
        # Dirty Equilibration using NVT and Alanine constraint
        #=============================================================================================

        self.simulation.context.setPositions(self.pdb_file.positions)
        
        system = self.simulation.system
        simulation = self.simulation
                
        nequib_steps = 5 #number of nvt equilibration steps with position constraints on Alanine
        Alanine_atoms = 22
                
        Alanine_masses = np.zeros(Alanine_atoms, np.double)
        Alanine_masses = Quantity(Alanine_masses, dalton)
        for i in range(Alanine_atoms):
            Alanine_masses[i]= system.getParticleMass(i)
            system.setParticleMass(i, 0.0)
        
        simulation.step(nequib_steps)
        
        for i in range(Alanine_atoms):
            system.setParticleMass(i,Alanine_masses[i].value_in_unit(dalton))
