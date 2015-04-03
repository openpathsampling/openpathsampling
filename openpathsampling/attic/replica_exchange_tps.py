'''
Created on 01.07.2014

@author: JD Chodera
@author ASJS Mey
@author: ...
'''

import copy
import numpy
import time
import os
import math

import netCDF4 as netcdf # for netcdf interface provided by netCDF4 in enthought python

from simtk.unit import nanosecond, picosecond, nanometers, nanometer, picoseconds, femtoseconds, femtosecond, kilojoules_per_mole, Quantity

from trajectory import Trajectory
from snapshot import Snapshot


#=============================================================================================
# SOURCE CONTROL
#=============================================================================================

__version__ = "$Id: KobAndersen.py 524 2010-01-05 07:47:29Z jchodera $"


#=============================================================================================
# REPLICA-EXCHANGE TRANSITION PATH SAMPLING
#=============================================================================================

class ReplicaExchangeTPS(object):
    """
    Replica-exchange transition path sampling among different values of s-field.

    """

    def __init__(self, system, ensembles, trajectories, ncfilename):
        """
        Initialize replica-exchange transition path sampling.

        ARGUMENTS

        system (simtk.chem.openState) - system object
        volume (list of TransitionPathSampling objects)
        trajectories - trajectory or list of trajectories to initialize TPS with
        ncfilename - name of NetCDF file to create or resume from

        TODO

        Implement a more general method where we have generateTrajectory(x0) report P(x0,x1,...,xT) and divide by P(x0) as appropriate.

        """
        # Store state information.
        self.system = system
        self.natoms = system.getNumParticles()

        # Store TransitionPathSampling objects
        self.ensembles = ensembles
        self.nstates = len(ensembles)

        # Distribute trajectories (and segment energies).
        self.trajectories = list()
        try:
            # Test if trajectories is an indexable set of Trajectory objects.
            shape = trajectories[0][0].coordinates.shape # size of first snapshot of first trajectory
            self.trajectories = [ copy.deepcopy(trajectories[index % len(trajectories)]) for index in range(self.nstates) ]
        except:
            try:
                # Test if trajectories is a single Trajectory object.
                shape = trajectories[0].coordinates.shape
                self.trajectories = [ copy.deepcopy(trajectories) for index in range(self.nstates) ]
            except:
                raise Exception("Don't know what to do with 'trajectories' object; doesn't seem to be a Trajectory object or list of Trajectories.")

        # Store number of frames per trajectory.
        self.n_frames = len(self.trajectories[0])

        # Set default options.
        self.number_of_iterations = 50
        self.title = 'Replica-exchange simulation'
        self.store_filename = ncfilename
        self.verbose = True
        
        # Flag as uninitialized.
        self._initialized = False

        return

    def _initialize(self):
        """
        Initialize the simulation, and bind to a storage file.

        """
        if self._initialized:
            print "Replica-exchange simulation has already been initialized."
            raise Exception("Already initialized.")

        print "Initializing..."

        # Allocate storage.
        self.replica_states     = numpy.zeros([self.nstates], numpy.int32) # replica_states[i] is the state that replica i is currently at
        for state_index in range(self.nstates):
            self.replica_states[state_index] = state_index
        
        self.log_P_kl           = numpy.zeros([self.nstates,self.nstates], numpy.float32) # log_P_kl[k,l] is log probability of replica k in state l
        self.swap_Pij_accepted  = numpy.zeros([self.nstates, self.nstates], numpy.float32) # swap_Pij_accepted[i,j] is fraction of swap attempts between states i and j that were accepted during last mixing phase
        self.activities         = [ None for replica_index in range(self.nstates) ] # activities[i] is activity of replica i
        self.path_hamiltonians  = [ None for replica_index in range(self.nstates) ] # activities[i] is activity of replica i        

        # Check if netcdf file extists.
        if os.path.exists(self.store_filename) and (os.path.getsize(self.store_filename) > 0):
            # Resume from NetCDF file.
            self._resume_from_netcdf()
        else:
            # Initialize current iteration counter.
            self.iteration = 0
            
            # Compute energies of all alchemical replicas
            self._compute_activities()
            
            # Initialize NetCDF file.
            self._initialize_netcdf()

            # Store initial state.
            self._write_iteration_netcdf()
  
        # Signal that the class has been initialized.
        self._initialized = True

        return

    def _finalize(self):
        """
        Do anything necessary to clean up.

        """

        self.ncfile.close()

        return
    def equilibrium_run(self):
        """
        This is just supposed to do an equilibrium simmulation of the system and export all the information in order to plot the s=0 distribution without using TPS
        """
        if not self._initialized:
            self.initialize()

    def run(self):
        """
        Run the replica-exchange simulation.

        Any parameter changes that were made between object creation and calling this method become locked in at this point,
        and the object will create and bind to the store file.

        """

        # Make sure we've initialized everything and bound to a storage file before we begin execution.
        if not self._initialized:
            self._initialize()

        # Main loop
        while (self.iteration < self.number_of_iterations):
            start_time = time.time()
            if self.verbose: print "\nIteration %d / %d" % (self.iteration+1, self.number_of_iterations)

            # Attempt replica swaps to sample from equilibrium permuation of states associated with replicas.
            self._mix_replicas()

            # Propagate replicas.
            self._propagate_replicas()

            # Compute energies of all replicas at all states.
            self._compute_activities()

            # Write to storage file.
            self._write_iteration_netcdf()
            
            # Increment iteration counter.
            self.iteration += 1
            end_time = time.time()
            if self.verbose: print "Iteration took %.3f s" % (end_time - start_time)

        # Clean up and close storage files.
        self._finalize()

        return

    def _propagate_replicas(self):
        """
        Propagate all replicas.
        
        NOTES
        
        self.ensembles contain only snapshot from this ensemble. The replicas are moving between these and self.replica_states
        contains the information necessary. But: there is no way to reconstruct the path a single replica takes!

        TODO

        * Create a parallel implementation when OpenMM supports this.

        """

        # Propagate all replicas.
        for replica_index in range(self.nstates):
            state_index = self.replica_states[replica_index]
            trajectory = self.trajectories[replica_index]
            trajectory = self.ensembles[state_index].sampleTrajectory(trajectory)
            self.trajectories[replica_index] = trajectory

        return

    def _compute_activities(self):
        """
        Compute activities K[x(t)] of all replicas.
        
        """

        print "Computing activities..."

        print "Computing trajectory probabilities..."

        # Compute path Hamiltonians for all replicas.
        for replica_index in range(self.nstates):
            trajectory = self.trajectories[replica_index]
            self.path_hamiltonians[replica_index] = self.ensembles[replica_index].pathHamiltonian(trajectory)

        # Compute activities for all replicas.
        for replica_index in range(self.nstates):
            trajectory = self.trajectories[replica_index]
            self.activities[replica_index] = self.ensembles[replica_index].computeActivity(trajectory)
            #print "replica %5d, x0 = %42s, K = %32s" % (replica_index, str(trajectory[0].coordinates[0,:]), str(self.activities[replica_index] / self.ensembles[replica_index].K_reduced_unit))

        # Compute log biasing probabilities for all replicas in all states.
        print self.log_P_kl # DEBUG
        for replica_index in range(self.nstates):
            K = self.activities[replica_index]
            H = self.path_hamiltonians[replica_index]
            for state_index in range(self.nstates):
                s = self.ensembles[state_index].s
                beta = self.ensembles[state_index].beta
                self.log_P_kl[replica_index,state_index] = - beta * H - s * K
        print self.log_P_kl # DEBUG
        
        if self.verbose:
            print "states = "
            for replica_index in range(self.nstates):
                print "%6d" % self.replica_states[replica_index],
            print ""
            print "path Hamiltonians = "
            for replica_index in range(self.nstates):
                print "%6.3f" % (self.path_hamiltonians[replica_index] / self.ensembles[replica_index].H_reduced_unit),
            print ""            
            print "activities = "
            for replica_index in range(self.nstates):
                print "%6.3f" % (self.activities[replica_index] / self.ensembles[replica_index].K_reduced_unit),
            print ""                
        return

    def _mix_replicas(self):
        """
        Attempt exchanges between all replicas to enhance mixing.

        NOTES

        This function might be slow in pure Python, so it may be necessary to re-code this in a compiled language.
        We certainly don't want this function to take a substantial fraction of the iteration time.
        
        """

        start_time = time.time()

        print "self.nstates = %d" % self.nstates

        # Determine number of swaps to attempt to ensure thorough mixing.
        # TODO: Replace this with analytical result computed to guarantee sufficient mixing.
        nswap_attempts = self.nstates**5 # number of swaps to attempt
        
        # Allocate storage to keep track of mixing.
        Nij_proposed = numpy.zeros([self.nstates,self.nstates], numpy.float32) # Nij_proposed[i][j] is the number of swaps proposed between states i and j, prior of 1
        Nij_accepted = numpy.zeros([self.nstates,self.nstates], numpy.float32) # Nij_proposed[i][j] is the number of swaps proposed between states i and j, prior of 1

        # Show log P
        if self.verbose:
            print "log_P[replica,state] ="
            print self.log_P_kl # DEBUG            
            print "%6s" % "",
            for jstate in range(self.nstates):
                print "%6d" % jstate,
            print ""
            for ireplica in range(self.nstates):
                print "%-6d" % ireplica,
                for jstate in range(self.nstates):
                    log_P = self.log_P_kl[ireplica,jstate]
                    print "%8.3f" % log_P,
                print ""

        # Attempt swaps to mix replicas.
        nswaps_accepted = 0
        for swap_attempt in range(nswap_attempts):
            # Choose replicas to attempt to swap.
            i = numpy.random.randint(self.nstates) # Choose replica i uniformly from set of replicas.
            j = numpy.random.randint(self.nstates) # Choose replica j uniformly from set of replicas.

            # Determine which states these resplicas correspond to.
            istate = self.replica_states[i] # state in replica slot i
            jstate = self.replica_states[j] # state in replica slot j

            # Compute log probability of swap.
            log_P_accept = (self.log_P_kl[i,jstate] + self.log_P_kl[j,istate]) - (self.log_P_kl[i,istate] + self.log_P_kl[j,jstate])

            #print "replica (%3d,%3d) states (%3d,%3d) energies (%8.1f,%8.1f) %8.1f -> (%8.1f,%8.1f) %8.1f : log_P_accept %8.1f" % (i,j,istate,jstate,self.u_kl[i,istate],self.u_kl[j,jstate],self.u_kl[i,istate]+self.u_kl[j,jstate],self.u_kl[i,jstate],self.u_kl[j,istate],self.u_kl[i,jstate]+self.u_kl[j,istate],log_P_accept)
            
            # Record that this move has been proposed.
            Nij_proposed[istate,jstate] += 0.5
            Nij_proposed[jstate,istate] += 0.5

            # Accept or reject.
            if (log_P_accept >= 0.0 or (numpy.random.rand() < math.exp(log_P_accept))):
                # Swap states in replica slots i and j.
                (self.replica_states[i], self.replica_states[j]) = (self.replica_states[j], self.replica_states[i])
                # Accumulate statistics
                Nij_accepted[istate,jstate] += 0.5
                Nij_accepted[jstate,istate] += 0.5
                nswaps_accepted += 1

        # Report statistics of acceptance.
        swap_fraction_accepted = float(nswaps_accepted) / float(nswap_attempts);
  
        # Estimate transition probabilities between all states.
        swap_Pij_accepted = numpy.zeros([self.nstates,self.nstates], numpy.float32)
        for istate in range(self.nstates):
            for jstate in range(self.nstates):
                if (Nij_proposed[istate,jstate] > 0.0):
                    swap_Pij_accepted[istate,jstate] = Nij_accepted[istate,jstate] / Nij_proposed[istate,jstate]
                else:
                    swap_Pij_accepted[istate,jstate] = 0.0
        self.swap_Pij_accepted = swap_Pij_accepted

        end_time = time.time()

        # Report on mixing.
        # TODO: Add this behind a verbose flag.
        if self.verbose:
            PRINT_CUTOFF = 0.001 # Cutoff for displaying fraction of accepted swaps.
            print "Fraction of accepted swaps between states:"
            print "%6s" % "",
            for jstate in range(self.nstates):
                print "%6d" % jstate,
            print ""
            for istate in range(self.nstates):
                print "%-6d" % istate,
                for jstate in range(self.nstates):
                    P = self.swap_Pij_accepted[istate,jstate]
                    if (P >= PRINT_CUTOFF):
                        print "%6.3f" % P,
                    else:
                        print "%6s" % "",
                print ""
            print "Mixing of replicas took %.3f s" % (end_time - start_time)

        return


    def _initialize_netcdf(self):
        """
        Initialize NetCDF file for storage.
        
        """    

        # Open NetCDF file for writing
        # ncfile = netcdf.NetCDFFile(self.store_filename, 'w') # for Scientific.IO.NetCDF
        ncfile = netcdf.Dataset(self.store_filename, 'w') # for netCDF4

        # Create dimensions.
        ncfile.createDimension('iteration', 0) # unlimited number of iterations
        ncfile.createDimension('replica', self.nstates) # number of replicas
        ncfile.createDimension('frame', self.n_frames) # number of frames per trajectory
        ncfile.createDimension('atom', self.natoms) # number of atoms in system
        ncfile.createDimension('spatial', 3) # number of spatial dimensions

        # Set global attributes.
        setattr(ncfile, 'tile', self.title)
        setattr(ncfile, 'application', 'KobAndersen')
        setattr(ncfile, 'program', 'KobAndersen.py')
        setattr(ncfile, 'programVersion', __version__)
        setattr(ncfile, 'Conventions', 'ReplicaExchangeTPS')
        setattr(ncfile, 'ConventionVersion', '0.1')

        ensemble = self.ensembles[0]
        setattr(ncfile, 'sKfactor', ensemble.N * ensemble.t_obs / ensemble.delta_t) # the quantity by which s * K must be multipled to obtain log weight
        setattr(ncfile, 'betaHfactor', 1.0) # the quantity by which beta * H must be multipled to obtain log weight
        
        # Create variables.
        ncvar_fields      = ncfile.createVariable('fields', 'f', ('replica',)) # s fields (in reduced units)
        ncvar_fields      = ncfile.createVariable('betas', 'f', ('replica',)) # inverse temperatures (in reduced units)

        ncvar_trajectory_coordinates = ncfile.createVariable('trajectory_coordinates', 'f', ('replica','frame','atom','spatial'))
        ncvar_trajectory_velocities  = ncfile.createVariable('trajectory_velocities',  'f', ('replica','frame','atom','spatial'))
        ncvar_trajectory_potential   = ncfile.createVariable('trajectory_potential',   'f', ('replica','frame'))
        ncvar_trajectory_kinetic     = ncfile.createVariable('trajectory_kinetic',     'f', ('replica','frame'))

        ncvar_states      = ncfile.createVariable('states', 'i', ('iteration','replica'))
        ncvar_activities  = ncfile.createVariable('activities', 'f', ('iteration','replica'))
        ncvar_path_hamiltonians = ncfile.createVariable('path_hamiltonians', 'f', ('iteration','replica'))        
        ncvar_log_probabilities = ncfile.createVariable('log_probabilities', 'f', ('iteration','replica','replica'))
        ncvar_mixing      = ncfile.createVariable('mixing', 'f', ('iteration','replica','replica'))
        
        # Define units for variables.
        setattr(ncvar_trajectory_coordinates, 'units', 'nm')
        setattr(ncvar_trajectory_velocities,  'units', 'nm/ps')
        setattr(ncvar_trajectory_potential,   'units', 'kJ/mol')
        setattr(ncvar_trajectory_kinetic,     'units', 'kJ/mol')
        
        setattr(ncvar_states,    'units', 'none')
        setattr(ncvar_mixing,    'units', 'none')
        # TODO: fields and activities

        # Set display formatting attributes.
        #setattr(ncvar_trajectories, 'C_format', r'%9.4f')
        #setattr(ncvar_activities, 'C_format', r'%5.3f')
        
        # Define long (human-readable) names for variables.
        setattr(ncvar_states,    "long_uid", "states[iteration][replica] is the state index (0..nstates-1) of replica 'replica' of iteration 'iteration'.")
        setattr(ncvar_mixing,    "long_uid", "mixing[iteration][i][j] is the fraction of proposed transitions between states i and j that were accepted during mixing using the coordinates from iteration 'iteration-1'.")
        # TODO: other variables 

        # Store fields.
        for state_index in range(self.nstates):
            s = self.ensembles[state_index].s
            s_reduced_unit = self.ensembles[state_index].s_reduced_unit
            ncfile.variables['fields'][state_index] = s / s_reduced_unit

        # Store inverse temperatures.
        for state_index in range(self.nstates):
            beta = self.ensembles[state_index].beta
            beta_reduced_unit = self.ensembles[state_index].beta_reduced_unit
            ncfile.variables['betas'][state_index] = beta / beta_reduced_unit

        # Force sync to disk to avoid data loss.
        ncfile.sync()

        # Store netcdf file handle.
        self.ncfile = ncfile
        
        return
    
    def _write_iteration_netcdf(self):
        """
        Write positions, states, and energies of current iteration to NetCDF file.
        
        """

        # DEBUG
        ensemble = self.ensembles[0]
        setattr(self.ncfile, 'sKfactor', ensemble.N * ensemble.t_obs / ensemble.delta_t) # the quantity by which s * K must be multipled to obtain log weight
        setattr(self.ncfile, 'betaHfactor', 1.0) # the quantity by which beta * H must be multipled to obtain log weight

        # Store trajectories.
        for replica_index in range(self.nstates):
            trajectory = self.trajectories[replica_index]
            for frame_index in range(self.n_frames):                
                self.ncfile.variables['trajectory_coordinates'][replica_index,frame_index,:,:] = (trajectory[frame_index].coordinates / nanometers).astype(numpy.float32)
                self.ncfile.variables['trajectory_velocities'][replica_index,frame_index,:,:] = (trajectory[frame_index].velocities / (nanometers / picoseconds)).astype(numpy.float32)
                self.ncfile.variables['trajectory_potential'][replica_index,frame_index] = trajectory[frame_index].potential_energy / kilojoules_per_mole                                
                self.ncfile.variables['trajectory_kinetic'][replica_index,frame_index] = trajectory[frame_index].kinetic_energy / kilojoules_per_mole

        # Store state information.
        self.ncfile.variables['states'][self.iteration,:] = self.replica_states[:]

        # Store activities.
        for replica_index in range(self.nstates):
            state_index = self.replica_states[replica_index]
            K = self.activities[replica_index]
            K_reduced_unit = self.ensembles[state_index].K_reduced_unit
            self.ncfile.variables['activities'][self.iteration,replica_index] = K / K_reduced_unit

        # Store path Hamiltonians.
        for replica_index in range(self.nstates):
            state_index = self.replica_states[replica_index]
            H = self.path_hamiltonians[replica_index]
            H_reduced_unit = self.ensembles[state_index].H_reduced_unit
            self.ncfile.variables['path_hamiltonians'][self.iteration,replica_index] = H / H_reduced_unit

        # Store log probabilities.
        print "writing log_probabilities..." # DEBUG
        print self.log_P_kl # DEBUG
        self.ncfile.variables['log_probabilities'][self.iteration,:,:] = self.log_P_kl[:,:]
        print self.log_P_kl # DEBUG
        
        # Store mixing statistics.
        self.ncfile.variables['mixing'][self.iteration,:,:] = self.swap_Pij_accepted[:,:]

        # Force sync to disk to avoid data loss.
        self.ncfile.sync()

        return

    def _resume_from_netcdf(self):
        """
        Resume execution by reading current positions and energies from a NetCDF file.
        
        """

        # Open NetCDF file for reading
        # ncfile = netcdf.NetCDFFile(self.store_filename, 'r') # for Scientific.IO.NetCDF
        ncfile = netcdf.Dataset(self.store_filename, 'r') # for netCDF4
        
        # TODO: Perform sanity check on file before resuming

        # Get current dimensions.
        self.iteration = ncfile.variables['activities'].shape[0] - 1
        self.nstates = ncfile.variables['activities'].shape[1]
        self.n_frames = ncfile.variables['trajectory_coordinates'].shape[1]
        self.natoms = ncfile.variables['trajectory_coordinates'].shape[2]

        print "iteration = %d, nstates = %d, natoms = %d" % (self.iteration, self.nstates, self.natoms)

        # Restore trajectories.
        self.trajectories = list()
        for replica_index in range(self.nstates):
            trajectory = Trajectory()
            for frame_index in range(self.n_frames):                
                x = ncfile.variables['trajectory_coordinates'][replica_index,frame_index,:,:].astype(numpy.float32).copy()
                coordinates = Quantity(x, nanometers)                
                v = ncfile.variables['trajectory_velocities'][replica_index,frame_index,:,:].astype(numpy.float32).copy()
                velocities = Quantity(v, nanometers / picoseconds)                
                V = ncfile.variables['trajectory_potential'][replica_index,frame_index]
                potential_energy = Quantity(V, kilojoules_per_mole)
                T = ncfile.variables['trajectory_kinetic'][replica_index,frame_index]
                kinetic_energy = Quantity(T, kilojoules_per_mole)
                snapshot = Snapshot(coordinates=coordinates, velocities=velocities, kinetic_energy=kinetic_energy, potential_energy=potential_energy)
                trajectory.forward(snapshot)
            self.trajectories.forward(trajectory)

        # Restore state information.
        self.replica_states = ncfile.variables['states'][self.iteration,:].copy()

        # Restore log probabilities.
        print "Reading log probabilities..." # DEBUG
        print self.log_P_kl # DEBUG
        self.log_P_kl = ncfile.variables['log_probabilities'][self.iteration,:,:]
        print self.log_P_kl # DEBUG
        
        # Restore activities
        for replica_index in range(self.nstates):
            state_index = self.replica_states[replica_index]
            K_reduced_unit = self.ensembles[state_index].K_reduced_unit
            K = ncfile.variables['activities'][self.iteration,replica_index]
            self.activities[replica_index] = K * K_reduced_unit

        # Restore path Hamiltonians
        for replica_index in range(self.nstates):
            state_index = self.replica_states[replica_index]
            H_reduced_unit = self.ensembles[state_index].H_reduced_unit
            H = ncfile.variables['path_hamiltonians'][self.iteration,replica_index]
            self.path_hamiltonians[replica_index] = H * H_reduced_unit

        # Close NetCDF file.
        ncfile.close()        
        
        # Reopen NetCDF file for appending, and maintain handle.
        #self.ncfile = netcdf.NetCDFFile(self.store_filename, 'a') # for Scientific.IO.NetCDF
        self.ncfile = netcdf.Dataset(self.store_filename, 'a') # for netCDF4
        
        # DEBUG: Set number of iterations to be a bit more than we've done.
        #self.number_of_iterations = self.iteration + 50
        
        return