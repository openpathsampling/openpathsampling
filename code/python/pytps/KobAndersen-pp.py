#!/usr/local/bin/env python

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Replica-exchange transition path sampling in the s-field on a Kob-Andersen system.

This version has been modified to run in parallel using Parallel Python:

http://www.parallelpython.com/

DESCRIPTION


REFERENCES

[1] Hedges LO, Jack RL, Garrahan JP, and Chandler D. Dynamic order-disorder in atomic models
of structural glass-formers. Science 323:1309, 2009.

[2] Minh DDL and Chodera JD. Optimal estimators and asymptotic variances for nonequilibrium
path-ensemble averages. JCP 131:134110, 2009.

COPYRIGHT

@author John D. Chodera <jchodera@gmail.com>

This source file is released under the GNU General Public License.

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 
You should have received a copy of the GNU General Public License along with
this program.  If not, see <http://www.gnu.org/licenses/>.

"""

#=============================================================================================
# GLOBAL IMPORTS
#=============================================================================================

import os
import os.path
import sys
import numpy
import math
import copy
import time

import pp # parallel python

import scipy.optimize # necessary for systems with finnicky SciPy installations

import pyopenmm

import simtk
import simtk.chem.openmm as openmm
import simtk.unit as units

import simtk.chem.openmm.extras.optimize as optimize

#import Scientific.IO.NetCDF as netcdf # for netcdf interface in Scientific
import netCDF4 as netcdf # for netcdf interface provided by netCDF4 in enthought python

from trajectory import Snapshot, Trajectory, TransitionPathSampling

#=============================================================================================
# SOURCE CONTROL
#=============================================================================================

__version__ = "$Id: KobAndersen.py 524 2010-01-05 07:47:29Z jchodera $"

#=============================================================================================
# GLOBAL CONSTANTS
#=============================================================================================

kB = units.BOLTZMANN_CONSTANT_kB * units.AVOGADRO_CONSTANT_NA

#=============================================================================================
# Kob-Andersen two-component mixture of Lennard-Jones particles.
#=============================================================================================

def KobAndersen(N=150, NA=None, A_fraction=0.8, principal_component_density=0.96, mm=None, mass=None, epsilon=None, sigma=None, softcore=False, alpha=0.5, lambda_=1.0):
    """
    Create a test system containing a Kob-Andersen two-component mixture.

    A soft-core Lennard-Jones potential is used if 'softcore' is set to True.

    OPTIONAL ARGUMENTS

    N (int) - total number of atoms (default: 150)
    A_fraction (float) - fraction of A component
    principal_component_density (float) - NA sigma^3 / V (default: 0.96)
    softcore (bool) - soft-core Lennard Jones (Eq. 4 of Shirts and Pande, JCP 122:134508, 2005) is used if True (default: False)
    lambda_ (float) - alchemical parameter, where 1.0 is fully interacting, 0.0 is non-interacting (default: 1.0)
    alpha (float) - soft-core parameter (default: 0.5)

    RETURNS

    system (System)
    coordinates (numpy array)
    epsilon (simtk.unit) - fundamental energy scale (change to argument?)
    sigma (simtk.unit) - fundamental length scale (change to argument?

    EXAMPLES

    Create a Kob-Andersen two-component mixture.

    >>> epsilon     = 119.8 * units.kelvin * units.BOLTZMANN_CONSTANT_kB * units.AVOGADRO_CONSTANT_NA # arbitrary reference energy    
    >>> [system, coordinates] = KobAndersen(epsilon=epsilon)

    Create softcore Kob-Andersen two-component mixture with alchemical perturbation.

    >>> [system, coordinates] = KobAndersen(epsilon=epsilon, softcore=True, lambda_=0.0)

    Test the energy

    >>> # Create a Context.
    >>> kB = units.BOLTZMANN_CONSTANT_kB
    >>> NA = units.AVOGADRO_CONSTANT_NA
    >>> temperature = 0.6 * epsilon / kB / NA
    >>> collision_rate = 90.0 / units.picosecond
    >>> timestep = 1.0 * units.femtosecond    
    >>> integrator = openmm.LangevinIntegrator(temperature, collision_rate, timestep)
    >>> context = openmm.Context(system, integrator)
    >>> # Set positions
    >>> context.setPositions(coordinates)
    >>> # Evaluate the potential energy.
    >>> state = context.getState(getEnergy=True)
    >>> reduced_potential = (state.getPotentialEnergy() / epsilon)
    >>> print reduced_potential

    Integrate dynamics

    >>> nsteps = 1000 # number of steps to integrate
    >>> integrator.step(nsteps)
    >>> # Retrieve configuration to make sure no coordinates are nan
    >>> state = context.getState(getPositions=True)
    >>> coordinates = state.getPositions(asNumpy=True)
    >>> if numpy.any(numpy.isnan(coordinates / units.nanometers)): raise Exception('some coordinates are nan after integration: %s' % str(coordinates))

    """
    
    # Choose OpenMM package.
    if mm is None:
        mm = simtk.chem.openmm

    # Set unit system based on Rowley, Nicholson, and Parsonage argon parameters.
    if mass    is None:  mass        = 39.948 * units.amu # arbitrary reference mass        
    if epsilon is None:  epsilon     = 119.8 * units.kelvin * units.BOLTZMANN_CONSTANT_kB * units.AVOGADRO_CONSTANT_NA # arbitrary reference energy    
    if sigma   is None:  sigma       = 0.3405 * units.nanometers # arbitrary reference lengthscale

    # Define LJ mixture parameters.
    epsilon_AA  = 1.0 * epsilon
    epsilon_AB  = 1.5 * epsilon
    epsilon_BB  = 0.5 * epsilon

    sigma_AA    = 1.0 * sigma
    sigma_AB    = 0.8 * sigma
    sigma_BB    = 0.88 * sigma

    # Determine number of atoms of each component
    if (NA is not None):
        # User has specified number of A components
        NB = N - NA
    else:
        # Compute number of A components
        NA = int(math.floor(A_fraction * N))
        NB = N - NA
        
    # Create system
    system = mm.System()

    # Compute total system volume.
    volume = NA * sigma**3 / principal_component_density
    
    # Make system cubic in dimension.
    length = volume**(1./3.)
    # TODO: Can we change this to use tuples or 3x3 array?
    a = units.Quantity(numpy.array([1.0, 0.0, 0.0], numpy.float32), units.nanometer) * length/units.nanometer
    b = units.Quantity(numpy.array([0.0, 1.0, 0.0], numpy.float32), units.nanometer) * length/units.nanometer
    c = units.Quantity(numpy.array([0.0, 0.0, 1.0], numpy.float32), units.nanometer) * length/units.nanometer
    print "box edge length = %s" % str(length)
    system.setDefaultPeriodicBoxVectors(a, b, c)

    # Add particles to system.
    for n in range(NA):
        system.addParticle(mass)
    for n in range(NB):
        system.addParticle(mass)
            
    # Create nonbonded force term implementing Kob-Andersen two-component Lennard-Jones interaction.
    energy_expression = ""
    if not softcore:
        # Standard Lennard-Jones, truncated to 2.5 sigma and shifted so potential is zero at cutoff.
        energy_expression += '4.0*epsilon*((sigma/r)^12 - (sigma/r)^6) * step(2.5*sigma - r) - 4.0*epsilon*((1.0/2.5)^12 - (1.0/2.5)^6);'
    else:
        # Soft-core Lennard-Jones from Eq. 4 of Shirts and Pande, JCP 122:134508, 2005.
        # U_LJ(r) = \lambda 4 \epsilon_{ij} ( [alpha (1-lambda) + (r/sigma_{ij})^6]^(-2) - [alpha (1-lambda) + (r/sigma_{ij}^6]^-1 )
        energy_expression += '4.0*epsilon*lambda*inv*(inv - 1.0);'
        energy_expression += 'inv = (alpha*(1.0-lambda) + (r/sigma)^6)^(-1);' 
        
    # Add mixing rules for two types.
    energy_expression += "epsilon = epsilon0*(1.0*AA + 1.5*AB + 0.5*BB);"
    energy_expression += "sigma = sigma0*(1.0*AA + 0.8*AB + 0.88*BB);"
    energy_expression += "AB = 1.0 - AA - BB;"
    energy_expression += "AA = A1*A2;"
    energy_expression += "BB = B1*B2;"

    force = mm.CustomNonbondedForce(energy_expression)

    # Add alchemical global parameters.
    if softcore:
        force.addGlobalParameter('alpha', alpha)
        force.addGlobalParameter('lambda', lambda_)

    # Set epsilon0 and sigma0 global parameters.
    force.addGlobalParameter('epsilon0', epsilon)
    force.addGlobalParameter('sigma0', sigma)

    # Add per-particle parameters to indicate whether each particle is type A or B.
    force.addPerParticleParameter('A')
    force.addPerParticleParameter('B')

    # Add A and B particle identifications.
    for n in range(NA):
        force.addParticle((1.0, 0.0))
    for n in range(NB):
        force.addParticle((0.0, 1.0))
    
    # Set periodic boundary conditions with cutoff.
    force.setNonbondedMethod(mm.CustomNonbondedForce.CutoffPeriodic)
    cutoff = 2.5 * sigma_AA
    cutoff = min(cutoff, length/2.0*0.9999)
    print "setting cutoff distance to %s" % str(cutoff)
    print "twice cutoff is %s" % str(cutoff*2)
    force.setCutoffDistance(cutoff)    

    # Add nonbonded force term to the system.
    system.addForce(force)

    # Create initial coordinates using a Sobol' subrandom sequence in three dimensions.
    coordinates = numpy.zeros([N,3], numpy.float32)
    for n in range(N):
        coordinates[n,:] = numpy.random.rand(3)
    coordinates = units.Quantity(coordinates, units.nanometer) * (length / units.nanometer)
       
    # Return system and coordinates.
    return (system, coordinates)

#=============================================================================================
# UTILITY FUNCTIONS
#=============================================================================================

def write_pdb(filename, trajectory, atoms):
    """
    Write out a trajectory as a multi-model PDB files.

    ARGUMENTS
       filename (string) - name of PDB file to be written
       trajectory (Trajectory)
       atoms (list of dict) - parsed PDB file ATOM entries from read_pdb() or construct_atom_list(); coordinates will be modified
    """

    # Create file.
    outfile = open(filename, 'w')

    nframes = len(trajectory)

    # Write trajectory as models
    for frame_index in range(nframes):
        outfile.write("MODEL     %4d\n" % (frame_index+1))
        coordinates = trajectory[frame_index].coordinates
        # Write ATOM records.
        for (index, atom) in enumerate(atoms):
            atom["x"] = "%8.3f" % (coordinates[index,0] / units.angstroms)
            atom["y"] = "%8.3f" % (coordinates[index,1] / units.angstroms)
            atom["z"] = "%8.3f" % (coordinates[index,2] / units.angstroms)
            outfile.write('ATOM  %(serial)5s %(atom)4s%(altLoc)c%(resName)3s %(chainID)c%(Seqno)5s   %(x)8s%(y)8s%(z)8s\n' % atom)

        outfile.write("ENDMDL\n")
        
    # Close file.
    outfile.close()

    return

def construct_atom_list(N, NA):
    """
    Construct atom metadata for a Kob-Andersen system containing NA atoms of type A and (N-NA) atoms of type B.

    ARGUMENTS
       N (int) - total number of atoms
       NA (int) - number of A type atoms

    RETURN VALUES
       atoms (list of dict) - list of atom dictionaries with PDB metadata

    """

    NB = N - NA
    atoms = list()

    index = 1
    for n in range(NA):
        atom = dict()
        atom['serial'] = index
        atom['atom'] = ' Ar '
        atom['altLoc'] = ' '
        atom['resName'] = 'Ar '
        atom['chainID'] = ' '
        atom['Seqno'] = '%5d' % index        
        index += 1
        atoms.append(atom)        
    for n in range(NB):
        atom = dict()
        atom['serial'] = index
        atom['atom'] = ' He '
        atom['altLoc'] = ' '
        atom['resName'] = 'He '
        atom['chainID'] = ' '
        atom['Seqno'] = '%5d' % index        
        index += 1
        atoms.append(atom)
        
    return atoms

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

        system (simtk.chem.openmm.State) - system object
        ensemble (list of TransitionPathSampling objects)
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
        self.nframes = len(self.trajectories[0])

        # Set default options.
        self.number_of_iterations = 50
        self.title = 'Replica-exchange simulation'
        self.store_filename = ncfilename
        self.verbose = True        
        self.job_server = None
        self.ncpus = None

        # Flag as uninitialized.
        self._initialized = False

        return
    
    def _initialize(self):
        """
        Initialize the simulation, and bind to a storage file.

        """
        if self._initialized:
            print "Replica-exchange simulation has already been initialized."
            raise Error

        print "Initializing..."

        # Set up parallel Python work pool.
#        if not self.job_server:
#            if self.verbose: print "Setting up parallel Python worker pool..."
#            ppservers = ("*",)
#            #self.job_server = pp.Server()
#            self.job_server = pp.Server(ncpus=2,ppservers=ppservers)            
#            if self.ncpus is not None:
#                self.job_server.set_ncpus(ncpus)
#            if self.verbose: print "Starting %d workers on server" % (self.job_server.get_ncpus())
#            if self.verbose: print self.job_server.get_active_nodes()

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

        TODO

        * Create a parallel implementation when OpenMM supports this.

        """

        if self.job_server is None:
            # Update trajectory samples with TPS serially.
            for replica_index in range(self.nstates):
                state_index = self.replica_states[replica_index]
                trajectory = self.trajectories[replica_index]
                trajectory = self.ensembles[state_index].sampleTrajectory(trajectory)
                self.trajectories[replica_index] = trajectory
        else:
            # Update trajectory samples with TPS in parallel.
            jobs = dict()
            if self.verbose: print "Submitting TPS jobs to job server..."
            for replica_index in range(self.nstates):
                state_index = self.replica_states[replica_index]
                ensemble = self.ensembles[state_index]
                trajectory = self.trajectories[replica_index]
                ngpus = 2 # DEBUG
                deviceid = replica_index % ngpus
                #job = self.job_server.submit(ensemble.sampleTrajectory, (trajectory,), (), ('UserList', 'trajectory', 'trajectory.Snapshot', 'trajectory.Trajectory', 'trajectory.TransitionPathSampling', 'simtk.unit'))
                job = self.job_server.submit(ensemble.sampleTrajectory, (trajectory,deviceid), (), ('trajectory',))    
                jobs[replica_index] = job

            if self.verbose: print "Waiting for jobs to finish..."
            for replica_index in range(self.nstates):
                job = jobs[replica_index] # get job handle
                trajectory = job() # retrieve new trajectory
                self.trajectories[replica_index] = trajectory 
            if self.verbose: print "Done...."
                
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
        ncfile.createDimension('frame', self.nframes) # number of frames per trajectory
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
        setattr(ncvar_states,    "long_name", "states[iteration][replica] is the state index (0..nstates-1) of replica 'replica' of iteration 'iteration'.")
        setattr(ncvar_mixing,    "long_name", "mixing[iteration][i][j] is the fraction of proposed transitions between states i and j that were accepted during mixing using the coordinates from iteration 'iteration-1'.")
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
            for frame_index in range(self.nframes):                
                self.ncfile.variables['trajectory_coordinates'][replica_index,frame_index,:,:] = (trajectory[frame_index].coordinates / units.nanometers).astype(numpy.float32)
                self.ncfile.variables['trajectory_velocities'][replica_index,frame_index,:,:] = (trajectory[frame_index].velocities / (units.nanometers / units.picoseconds)).astype(numpy.float32)
                self.ncfile.variables['trajectory_potential'][replica_index,frame_index] = trajectory[frame_index].potential_energy / units.kilojoules_per_mole                                
                self.ncfile.variables['trajectory_kinetic'][replica_index,frame_index] = trajectory[frame_index].kinetic_energy / units.kilojoules_per_mole

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
        self.nframes = ncfile.variables['trajectory_coordinates'].shape[1]
        self.natoms = ncfile.variables['trajectory_coordinates'].shape[2]

        print "iteration = %d, nstates = %d, natoms = %d" % (self.iteration, self.nstates, self.natoms)

        # Restore trajectories.
        self.trajectories = list()
        for replica_index in range(self.nstates):
            trajectory = Trajectory()
            for frame_index in range(self.nframes):                
                x = ncfile.variables['trajectory_coordinates'][replica_index,frame_index,:,:].astype(numpy.float32).copy()
                coordinates = units.Quantity(x, units.nanometers)                
                v = ncfile.variables['trajectory_velocities'][replica_index,frame_index,:,:].astype(numpy.float32).copy()
                velocities = units.Quantity(v, units.nanometers / units.picoseconds)                
                V = ncfile.variables['trajectory_potential'][replica_index,frame_index]
                potential_energy = units.Quantity(V, units.kilojoules_per_mole)
                T = ncfile.variables['trajectory_kinetic'][replica_index,frame_index]
                kinetic_energy = units.Quantity(T, units.kilojoules_per_mole)
                snapshot = Snapshot(coordinates=coordinates, velocities=velocities, kinetic_energy=kinetic_energy, potential_energy=potential_energy)
                trajectory.append(snapshot)
            self.trajectories.append(trajectory)

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

#=============================================================================================
# MAIN AND TESTS
#=============================================================================================

def driver():
    # Constant
    kB = units.BOLTZMANN_CONSTANT_kB * units.AVOGADRO_CONSTANT_NA

    # Create a Kob-Andersen two-component mixture.
    N = 150 # number of particles
    A_fraction = 0.8 # fraction of A component
    NA = int(math.floor(A_fraction * N)) # number of A component
    mass        = 39.948 * units.amu # arbitrary reference mass        
    epsilon     = 119.8 * units.kelvin * units.BOLTZMANN_CONSTANT_kB * units.AVOGADRO_CONSTANT_NA # arbitrary reference energy    
    sigma       = 0.3405 * units.nanometers # arbitrary reference lengthscale    
    [system, coordinates] = KobAndersen(N=N, NA=NA, mass=mass, epsilon=epsilon, sigma=sigma)    

    # Choose the platform and device we will run on.
    #platform = openmm.Platform.getPlatformByName('Reference') # slow CPU platform
    platform_name = 'OpenCL'
    platform = openmm.Platform.getPlatformByName(platform_name) # fast GPU platform
    platform.setPropertyDefaultValue('OpenCLDeviceIndex', '0') # use first OpenCL device

    # Select test or production mode.
    # Test mode runs more quickly, but runs with a different set of parameters.
    #mode = 'test'
    mode = 'production'
    
    # Relevant times and timescales
    reduced_time = units.sqrt((mass * sigma**2) / (48.0 * epsilon)) # factor on which all reduced times are based
    timestep = 0.035 * reduced_time # velocity Verlet timestep from Ref. [1], Supplementary Section 1.1
    delta_t = (40.0 / 3.0) * reduced_time # number of timesteps per Delta t in Lester's paper, Ref. [1] Supplementary Section 1.1, specified exactly by Lester in private communication

    if mode == 'test':
        quenched_temperature = 0.6 * epsilon / kB # temperature corresponding to left column of Fig. 2 from [1]
        tau = 200.0 * reduced_time # private communication from Lester Hedges (tau = 60 * reduced time for Tred = 0.7, 200 * reduced_time for Tred = 0.6)
        t_obs = 1 * tau # number of intervals of tau, from Fig. 2 of [1] (red dashed line, right column,) # should be 33 * tau for Tred = 0.7, 20*tau for Tred = 0.6
    elif mode == 'production':
        quenched_temperature = 0.7 * epsilon / kB # temperature corresponding to left column of Fig. 2 from [1]
        tau = 60.0 * reduced_time # private communication from Lester Hedges (tau = 60 * reduced time for Tred = 0.7, 200 * reduced_time for Tred = 0.6)
        t_obs = 33 * tau # number of intervals of tau, from Fig. 2 of [1] (red dashed line, right column,) # should be 33 * tau for Tred = 0.7, 20*tau for Tred = 0.6        
    else:
        raise Exception("Unknown mode: %s" % mode)

    print "UNCORRECTED TIMES"
    print "timestep = %s" % str(timestep.in_units_of(units.femtosecond))
    print "delta_t = %s (%f steps)" % (str(delta_t.in_units_of(units.picosecond)), delta_t / timestep)
    print "tau = %s (%f steps)" % (str(tau.in_units_of(units.picosecond)), tau / timestep)
    print "t_obs = %s (%f steps)" % (str(t_obs.in_units_of(units.picosecond)), t_obs / timestep)

    # Determine integral number of steps per velocity randomization and number of trajectories
    nsteps_per_frame = int(round(delta_t / timestep)) # number of steps per trajectory frame and velocity randomization
    nframes = int(round(t_obs / delta_t)) # number of frames per trajectory
    print "number of steps per trajectory frame and velocity randomization = %d" % nsteps_per_frame
    print "number of frames per trajectory = %d" % nframes
    # Correct delta_t and t_obs to be integral
    delta_t = nsteps_per_frame * timestep 
    t_obs = nframes * delta_t
    print "CORRECTED TIMES (after rounding number of steps to integers)"
    print "timestep = %s" % str(timestep.in_units_of(units.femtosecond))
    print "delta_t = %s (%f steps)" % (str(delta_t.in_units_of(units.picosecond)), delta_t / timestep)
    print "tau = %s (%f steps)" % (str(tau.in_units_of(units.picosecond)), tau / timestep)
    print "t_obs = %s (%f steps)" % (str(t_obs.in_units_of(units.picosecond)), t_obs / timestep)

    # Replica-exchange filename
    ncfilename = 'repex.nc'

    # Set default simulation protocol.
    minimize = True
    equilibrate = True
    quench = True
    seed = True
    
    # If we're resuming from a previously-existing NetCDF file, change the protocol to skip equilibration and quenching.
    if os.path.exists(ncfilename):
        # No need to do these things if we are resuming
        minimize = True
        equilibrate = False
        quench = False
        seed = True
                
    if minimize:
        # Minimize the system prior to dynamics.
        print "Minimizing with L-BFGS..."
        # Initialize a minimizer with default options.
        minimizer = optimize.LBFGSMinimizer(system, verbose=True, platform=platform)
        # Minimize the initial coordinates.
        coordinates = minimizer.minimize(coordinates)
        # Clean up to release the Context.
        del minimizer

    # Set temperature for equilibration simulations.
    elevated_temperature = 2.0 * epsilon / kB

    if equilibrate:
        # Equilibrate at high temperature.
        collision_rate = 1.0 / delta_t
        nsteps = int(math.floor(t_obs / timestep))
        print "Equilibrating at %s for %d steps..." % (str(elevated_temperature), nsteps)
        integrator = openmm.LangevinIntegrator(elevated_temperature, collision_rate, timestep)
        context = openmm.Context(system, integrator, platform)
        context.setPositions(coordinates)
        integrator.step(nsteps)
        state = context.getState(getPositions=True)
        coordinates = state.getPositions(asNumpy=True)
        del state, context, integrator

    if quench:
        # Quench to final temperature
        collision_rate = 1.0 / delta_t
        nsteps = int(math.floor(t_obs / timestep))
        print "Quenching to %s for %d steps..." % (str(quenched_temperature), nsteps)
        integrator = openmm.LangevinIntegrator(quenched_temperature, collision_rate, timestep)
        context = openmm.Context(system, integrator, platform)
        context.setPositions(coordinates)
        integrator.step(nsteps)
        state = context.getState(getPositions=True)
        coordinates = state.getPositions(asNumpy=True)
        del state, context, integrator

    # Create a number of transition path sampling ensembles at different values of the field parameter s.
    s_reduced_unit = 1.0 / (sigma**2 * delta_t)
    #svalues = [0.00, 0.01, 0.02, 0.03, 0.04, 0.06] # minimal set of s values
    svalues = [0.00, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, -0.005]
    print "Initializing transition path sampling ensembles..."
    ensembles = [ TransitionPathSampling(system, N, NA, mass, epsilon, sigma, timestep, nsteps_per_frame, nframes, quenched_temperature, s * s_reduced_unit, platform_name=platform_name) for s in svalues ]

    trajectory = None
    if seed:
        # Generate an initial trajectory at zero field
        print "Generating seed trajectory for TPS..."
        trajectory = ensembles[0].generateTrajectory(coordinates, nframes)

    # Initialize replica-exchange TPS simulation.
    print "Initializing replica-exchange TPS..."
    simulation = ReplicaExchangeTPS(system, ensembles, trajectory, ncfilename)

    # Set simulation parameters.
    simulation.number_of_iterations = 2000

    # Run simulation
    print "Running replica-exchange TPS..."
    simulation.run()
        
    print "Work Done :) Go home and have a beer."

#=============================================================================================
# This driver sets up a simulation for multiple temperatures and s-values.
#=============================================================================================
def multitemp_driver():
    # Constant
    kB = units.BOLTZMANN_CONSTANT_kB * units.AVOGADRO_CONSTANT_NA

    # Create a Kob-Andersen two-component mixture.
    N = 150 # number of particles
    A_fraction = 0.8 # fraction of A component
    NA = int(math.floor(A_fraction * N)) # number of A component
    mass        = 39.948 * units.amu # arbitrary reference mass        
    epsilon     = 119.8 * units.kelvin * units.BOLTZMANN_CONSTANT_kB * units.AVOGADRO_CONSTANT_NA # arbitrary reference energy    
    sigma       = 0.3405 * units.nanometers # arbitrary reference lengthscale    
    [system, coordinates] = KobAndersen(N=N, NA=NA, mass=mass, epsilon=epsilon, sigma=sigma)    

    # Choose the platform and device we will run on.
    #platform = openmm.Platform.getPlatformByName('Reference') # slow CPU platform
    platform_name = 'OpenCL'
    platform = openmm.Platform.getPlatformByName(platform_name) # fast GPU platform
    platform.setPropertyDefaultValue('OpenCLDeviceIndex', '0') # use first OpenCL device

    # Select test or production mode.
    # Test mode runs more quickly, but runs with a different set of parameters.
    #mode = 'test'
    mode = 'production'
    
    # Relevant times and timescales
    reduced_time = units.sqrt((mass * sigma**2) / (48.0 * epsilon)) # factor on which all reduced times are based
    timestep = 0.035 * reduced_time # velocity Verlet timestep from Ref. [1], Supplementary Section 1.1
    delta_t = (40.0 / 3.0) * reduced_time # number of timesteps per Delta t in Lester's paper, Ref. [1] Supplementary Section 1.1, specified exactly by Lester in private communication

    if mode == 'test':
        quenched_temperature = 0.6 * epsilon / kB # temperature corresponding to left column of Fig. 2 from [1]
        tau = 200.0 * reduced_time # private communication from Lester Hedges (tau = 60 * reduced time for Tred = 0.7, 200 * reduced_time for Tred = 0.6)
        t_obs = 1 * tau # number of intervals of tau, from Fig. 2 of [1] (red dashed line, right column,) # should be 33 * tau for Tred = 0.7, 20*tau for Tred = 0.6
    elif mode == 'production':
        quenched_temperature = 0.7 * epsilon / kB # temperature corresponding to left column of Fig. 2 from [1]
        tau = 60.0 * reduced_time # private communication from Lester Hedges (tau = 60 * reduced time for Tred = 0.7, 200 * reduced_time for Tred = 0.6)
        t_obs = 33 * tau # number of intervals of tau, from Fig. 2 of [1] (red dashed line, right column,) # should be 33 * tau for Tred = 0.7, 20*tau for Tred = 0.6        
    else:
        raise Exception("Unknown mode: %s" % mode)

    print "UNCORRECTED TIMES"
    print "timestep = %s" % str(timestep.in_units_of(units.femtosecond))
    print "delta_t = %s (%f steps)" % (str(delta_t.in_units_of(units.picosecond)), delta_t / timestep)
    print "tau = %s (%f steps)" % (str(tau.in_units_of(units.picosecond)), tau / timestep)
    print "t_obs = %s (%f steps)" % (str(t_obs.in_units_of(units.picosecond)), t_obs / timestep)

    # Determine integral number of steps per velocity randomization and number of trajectories
    nsteps_per_frame = int(round(delta_t / timestep)) # number of steps per trajectory frame and velocity randomization
    nframes = int(round(t_obs / delta_t)) # number of frames per trajectory
    print "number of steps per trajectory frame and velocity randomization = %d" % nsteps_per_frame
    print "number of frames per trajectory = %d" % nframes
    # Correct delta_t and t_obs to be integral
    delta_t = nsteps_per_frame * timestep 
    t_obs = nframes * delta_t
    print "CORRECTED TIMES (after rounding number of steps to integers)"
    print "timestep = %s" % str(timestep.in_units_of(units.femtosecond))
    print "delta_t = %s (%f steps)" % (str(delta_t.in_units_of(units.picosecond)), delta_t / timestep)
    print "tau = %s (%f steps)" % (str(tau.in_units_of(units.picosecond)), tau / timestep)
    print "t_obs = %s (%f steps)" % (str(t_obs.in_units_of(units.picosecond)), t_obs / timestep)

    # Replica-exchange filename
    ncfilename = 'repex-Ts.nc' # DEBUG

    # Set default simulation protocol.
    minimize = True
    equilibrate = True
    quench = True
    seed = True
    
    # If we're resuming from a previously-existing NetCDF file, change the protocol to skip equilibration and quenching.
    if os.path.exists(ncfilename):
        # No need to do these things if we are resuming
        minimize = False
        equilibrate = False
        quench = False
        seed = False
                
    if minimize:
        # Minimize the system prior to dynamics.
        print "Minimizing with L-BFGS..."
        # Initialize a minimizer with default options.
        minimizer = optimize.LBFGSMinimizer(system, verbose=True, platform=platform)
        # Minimize the initial coordinates.
        coordinates = minimizer.minimize(coordinates)
        # Clean up to release the Context.
        del minimizer

    # Set temperature for equilibration simulations.
    elevated_temperature = 2.0 * epsilon / kB

    if equilibrate:
        # Equilibrate at high temperature.
        collision_rate = 1.0 / delta_t
        nsteps = int(math.floor(t_obs / timestep))
        print "Equilibrating at %s for %d steps..." % (str(elevated_temperature), nsteps)
        integrator = openmm.LangevinIntegrator(elevated_temperature, collision_rate, timestep)
        context = openmm.Context(system, integrator, platform)
        context.setPositions(coordinates)
        integrator.step(nsteps)
        state = context.getState(getPositions=True)
        coordinates = state.getPositions(asNumpy=True)
        del state, context, integrator

    if quench:
        # Quench to final temperature
        collision_rate = 1.0 / delta_t
        nsteps = int(math.floor(t_obs / timestep))
        print "Quenching to %s for %d steps..." % (str(quenched_temperature), nsteps)
        integrator = openmm.LangevinIntegrator(quenched_temperature, collision_rate, timestep)
        context = openmm.Context(system, integrator, platform)
        context.setPositions(coordinates)
        integrator.step(nsteps)
        state = context.getState(getPositions=True)
        coordinates = state.getPositions(asNumpy=True)
        del state, context, integrator

    # Create a number of transition path sampling ensembles at different values of the field parameter s and different temperatures
    s_reduced_unit = 1.0 / (sigma**2 * delta_t) # reduced units for s field
    T_reduced_unit = epsilon / kB # reduced units for temperature

    svalues = [0.00, 0.01, 0.02, 0.03, 0.04, 0.06] # minimal set of s values
    Tvalues = [0.70, 0.80, 0.90, 1.00, 1.10, 1.20] # some made-up temperatures

    # Create a pure Python representation of the openmm Swig System object for transport over the network.
    pyopenmm_system = pyopenmm.System(system)
               
    print "Initializing transition path sampling ensembles..."
    ensembles = [ TransitionPathSampling(pyopenmm_system, N, NA, mass, epsilon, sigma, timestep, nsteps_per_frame, nframes, T * T_reduced_unit, s * s_reduced_unit, platform_name=platform_name) for (s,T) in zip(svalues,Tvalues) ]

    trajectory = None
    if seed:
        # Generate an initial trajectory at zero field
        print "Generating seed trajectory for TPS..."
        trajectory = ensembles[0].generateTrajectory(coordinates, nframes)
    else:
        trajectory = list()
        for frame in range(nframes):
            trajectory.append(Snapshot(coordinates=coordinates))

    # Initialize replica-exchange TPS simulation.
    print "Initializing replica-exchange TPS..."
    simulation = ReplicaExchangeTPS(pyopenmm_system, ensembles, trajectory, ncfilename)

    # Set simulation parameters.
    simulation.number_of_iterations = 10

    # Run simulation
    print "Running replica-exchange TPS..."
    simulation.run()

    # Done.
    print "Work done!  Go home and mix yourself a cocktail."

if __name__ == "__main__":
    #import doctest
    #doctest.testmod()    

    # Run replica-exchange TPS driver
    #driver()

    # Run new-and-improved multiple temperature and s value driver.
    multitemp_driver()

    
   
  
