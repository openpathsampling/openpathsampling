#!/usr/local/bin/env python

print "importing mpi4py.MPI..."
from mpi4py import MPI
print "Started node %d / %d" % (MPI.COMM_WORLD.rank, MPI.COMM_WORLD.size)

#=============================================================================================
# MODULE DOCSTRING
#=============================================================================================

"""
Replica-exchange transition path sampling in the s-field on a Kob-Andersen system.
Parallelization using mpi4py.

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

import scipy.optimize # necessary for systems with finnicky SciPy installations

import simtk
import simtk.chem.openmm as openmm
import simtk.unit as units

import simtk.chem.openmm.extras.optimize as optimize

#import Scientific.IO.NetCDF as netcdf # for netcdf interface in Scientific
import netCDF4 as netcdf # for netcdf interface provided by netCDF4 in enthought python

#from mpi4py import MPI # MPI wrapper

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
# SIMULATION SNAPSHOT 
#=============================================================================================

class Snapshot(object):
    """
    Simulation snapshot.

    """
    #def __init__(self, context=None, coordinates=None, velocities=None, box_vectors=None, potential_energy=None, kinetic_energy=None):
    def __init__(self, context=None, coordinates=None, box_vectors=None, potential_energy=None, kinetic_energy=None):        
        """
        Create a simulation snapshot from either an OpenMM context or individually-specified components.

        OPTIONAL ARGUMENTS

        context (simtk.chem.openmm.Context) - if not None, the current state will be queried to populate simulation snapshot; otherwise, can specify individual components (default: None)
        coordinates (simtk.unit.Quantity wrapping Nx3 numpy array of dimension length) - atomic coordinates (default: None)
        velocities (simtk.unit.Quantity wrapping Nx3 numpy array of dimension length) - atomic velocities (default: None)
        box_vectors - periodic box vectors (default: None)
        potential_energy (simtk.unit.Quantity of units energy/mole) - potential energy at current timestep (default: None)
        kinetic_energy (simtk.unit.Quantity of units energy/mole) - kinetic energy at current timestep (default: None)
        
        """

        if context is not None:
            # Get current state from OpenMM Context object.
            state = context.getState(getPositions=True, getVelocities=True, getEnergy=True)
            
            # Populate current snapshot data.
            self.coordinates = state.getPositions(asNumpy=True)
            #self.velocities = state.getVelocities(asNumpy=True)
            self.box_vectors = state.getPeriodicBoxVectors() # TODO: set asNumpy=True once bug in OpenMM is fixed
            self.potential_energy = state.getPotentialEnergy()
            self.kinetic_energy = state.getKineticEnergy()
        else:
            if coordinates is not None: self.coordinates = copy.deepcopy(coordinates)
            #if velocities is not None: self.velocities = copy.deepcopy(velocities)
            if box_vectors is not None: self.box_vectors = copy.deepcopy(box_vectors)
            if potential_energy is not None: self.potential_energy = copy.deepcopy(potential_energy)
            if kinetic_energy is not None: self.kinetic_energy = copy.deepcopy(kinetic_energy)                       

        # Check for nans in coordinates, and raise an exception if something is wrong.
        if numpy.any(numpy.isnan(self.coordinates)):
            raise Exception("Some coordinates became 'nan'; simulation is unstable or buggy.")

        return

    @property
    def total_energy(self):
        return self.kinetic_energy + self.potential_energy

#=============================================================================================
# SIMULATION TRAJECTORY
#=============================================================================================

class Trajectory(list):
    """
    Simulation trajectory.

    """
    def __init__(self, trajectory=None):
        """
        Create a simulation trajectory object

        OPTIONAL ARGUMENTS

        trajectory (Trajectory) - if specfiied, make a deep copy of specified trajectory
        
        """

        # Initialize list.
        list.__init__(self)

        if trajectory is not None:
            # Try to make a copy out of whatever container we were provided
            for snapshot in trajectory:
                snapshot_copy = copy.deepcopy(snapshot)                    
                self.append(snapshot_copy)

        return

    def reverse(self):
        """
        Reverse the trajectory.

        NOTE

        We cannot handle the velocities correctly when reversing the trajectory, so velocities will no longer be meaningful.
        Kinetic energies are correctly updated, however, and path actions should be accurate.

        """
        # Reverse the order of snapshots within the trajectory.
        list.reverse(self)

        # Determine number of snapshots.
        nsnapshots = self.__len__()
        
        # Recalculate kinetic energies for the *beginning* of each trajectory segment.
        # This makes use of the fact that the energy is (approximately) conserved over each trajectory segment, in between velocity randomizations.
        # Note that this may be a poor approximation in some cases.
        for t in range(nsnapshots-1):
            self[t].kinetic_energy = self[t+1].total_energy - self[t].potential_energy

        # No use reversing momenta, since we can't determine what appropriate reversed momenta should be.
        
        return
    
#=============================================================================================
# TRANSITION PATH SAMPLING
#=============================================================================================

class TransitionPathSampling(object):
    """
    Transition path sampling with fixed s-field.

    """

    def __init__(self, system, N, NA, mass, epsilon, sigma, timestep, nsteps_per_frame, nframes, temperature, s, platform=None, verbose=False):
        """
        Initialize a transition path sampling simulation with activity field strength s.

        ARGUMENTS

        system (simtk.chem.openmm.System) - the molecular mechanics system
        N (int) - number of particles in the system
        NA (int) - number of A-type particles to be used in computing activity functional K[x(t)]
        mass (simtk.unit.Quantity with units mass/particle) - particle mass to be used in computing reduced quantites 
        epsilon (simtk.unit.Quantity with units energy/mole) - Lennard-Jones well depth parameter to be used in computing reduced quantities
        sigma (simtk.unit.Quantity with units distance) - Lennard-Jones sigma parameter to be used in computing reduced quantities
        timestep (simtk.unit.Quantity with units time) - timestep to be used in dynamics simulations
        nsteps_per_frame (int) - number of MD timesteps per simulation snapshot or 'frame' in a trajectory
        nframes (int) - number of simulation snapshots or 'frames' for a fixed-length TPS trajectory, not counting initial configuration
        temperture (simtk.unit.Quantity with units temperature) - simulation temeprature
        s (float) - the value of the field parameter

        OPTIONAL ARGUMENTS

        platform (simtk.chem.openmm.Platform) - platform to use for OpenMM calculations
        
        """

        # Store local copy of System.
        self.system = system

        self.epsilon = epsilon
        self.sigma = sigma

        self.timestep = timestep
        self.nsteps_per_frame = nsteps_per_frame
        self.delta_t = timestep * nsteps_per_frame

        self.temperature = temperature

        # Compute thermal energy and inverse temperature from specified temperature.
        self.kT = kB * self.temperature # thermal energy
        self.beta = 1.0 / self.kT # inverse temperature

        # Store number of A particles.
        self.N = N
        self.NA = NA

        # Store field.
        self.s = s

        # Statistics
        self.nattempted = 0
        self.naccepted = 0

        # Create a Verlet integrator.
        self.integrator = openmm.VerletIntegrator(self.timestep)

        # Create a Context for integration.
        if platform:
            self.context = openmm.Context(self.system, self.integrator, platform)
        else:
            self.context = openmm.Context(self.system, self.integrator)            

        # Store reduced units
        self.t_obs = nframes * self.delta_t
        self.s_reduced_unit = 1.0 /  (self.sigma**2 * self.delta_t)
        self.K_reduced_unit = (self.N * self.t_obs * self.sigma**2)
        self.H_reduced_unit = self.epsilon # reduced unit for path Hamiltonian (energy)
        self.beta_reduced_unit = 1.0 / self.epsilon # reduced unit for inverse temperature

        # Form vectors of masses and sqrt(kT/m) for force propagation and velocity randomization.
        nparticles = system.getNumParticles()
        self.mass = units.Quantity(numpy.zeros([nparticles,3], numpy.float64), units.amu)
        for particle_index in range(nparticles):
            self.mass[particle_index,:] = self.system.getParticleMass(particle_index)
        kT = kB * temperature # thermal energy    
        self.sqrt_kT_over_m = units.Quantity(numpy.zeros([nparticles,3], numpy.float64), units.nanometers / units.picosecond)
        for particle_index in range(nparticles):
            self.sqrt_kT_over_m[particle_index,:] = units.sqrt(kT / self.mass[particle_index,0]) # standard deviation of velocity distribution for each coordinate for this atom

        self.verbose = verbose

        return

    def assignMaxwellBoltzmannVelocities(self, remove_com_velocity=False):
        """
        Generate Maxwell-Boltzmann velocities at the current simulation temperature.

        OPTIONAL ARGUMENTS

        remove_com_velocity (boolean) - if True, the center-of-mass velocity will be removed after the velocities are randomized (default: False)

        TODO

        This could be sped up by introducing vector operations.

        """

        # Get number of atoms
        nparticles = self.system.getNumParticles()

        # Assign velocities from the Maxwell-Boltzmann distribution.
        velocities = self.sqrt_kT_over_m * numpy.random.standard_normal(size=(nparticles,3))

        if remove_com_velocity:
            # Remove center of mass velocity
            velocity_units = self.sqrt_kT_over_m.unit
            com_velocity = units.Quantity(numpy.reshape((velocities / velocity_units).mean(0), (1,3)), velocity_units)
            velocities -= units.Quantity(numpy.repeat(com_velocity / velocity_units, nparticles, axis=0),velocity_units)

        # Return velocities
        return velocities

    def generateTrajectory(self, x0, nframes):
        """
        Generate a velocity Verlet trajectory consisting of ntau segments of tau_steps in between storage of Snapshots and randomization of velocities.

        ARGUMENTS
        
        x0 (coordinate set) - initial coordinates; velocities will be assigned from Maxwell-Boltzmann distribution
        nframes (int) - number of trajectory segments to generate

        RETURNS

        trajectory (list of Snapshot) - generated trajectory of initial conditions, including initial coordinate set

        NOTES

        This routine generates a velocity Verlet trajectory for systems without constraints by wrapping the OpenMM 'VerletIntegrator' in two half-kicks of the velocity.
        
        """

        # Set initial positions
        self.context.setPositions(x0)

        # Store initial state for each trajectory segment in trajectory.
        trajectory = Trajectory()

        # Generate trajectory segments.
        for frame_index in range(nframes):
            # Assign velocities from Maxwell-Boltzmann distribution
            velocities = self.assignMaxwellBoltzmannVelocities(remove_com_velocity=True)
            self.context.setVelocities(velocities)            

            # Store initial snapshot of trajectory segment.
            snapshot = Snapshot(context=self.context)
            trajectory.append(snapshot)
            
            # Propagate dynamics by velocity Verlet.
            # We only have leapfrog integrator available, so we wrap it in two half-kicks.
            # Back-kick by half a timestep to get ready for leapfrog integration.
            state = self.context.getState(getForces=True, getVelocities=True)
            force = state.getForces(asNumpy=True)
            velocities = state.getVelocities(asNumpy=True)
            velocities -= 0.5 * force/self.mass * self.timestep
            self.context.setVelocities(velocities)
            # Step using leapfrog.
            self.integrator.step(self.nsteps_per_frame)
            # Forward-kick by half a timestep to bring velocities into sync with positions.
            state = self.context.getState(getForces=True, getVelocities=True)
            force = state.getForces(asNumpy=True)
            velocities = state.getVelocities(asNumpy=True)
            velocities += 0.5 * force/self.mass * self.timestep            
            self.context.setVelocities(velocities)

        # Store final snapshot of trajectory.
        snapshot = Snapshot(self.context)
        trajectory.append(snapshot)

        return trajectory

    def logEquilibriumTrajectoryProbability(self, trajectory):
        """
        Compute the log equilibrium probability (up to an unknown additive constant) of an unbiased trajectory evolved according to Verlet dynamics with Andersen thermostatting.

        ARGUMENTS

        trajectory (Trajectory) - the trajectory

        RETURNS

        log_q (float) - the log equilibrium probability of the trajectory

        """

        nsnapshots = len(trajectory)
        log_q = - self.beta * trajectory[0].total_energy
        for snapshot_index in range(1, nsnapshots-1):
            log_q += - self.beta * trajectory[snapshot_index].kinetic_energy

        return log_q

    def pathHamiltonian(self, trajectory):
        """
        Compute the generalized path Hamiltonian of the trajectory.

        ARGUMENTS

        trajectory (Trajectory) - the trajectory

        RETURNS

        H (simtk.unit.Quantity with units of energy) - the generalized path Hamiltonian

        REFERENCES

        For a description of the path Hamiltonian, see [1]:

        [1] Chodera JD, Swope WC, Noe F, Prinz JH, Shirts MR, and Pande VS. Dynamical reweighting:
        Improved estimates of dynamical properties from simulations at multiple temperatures.    

        """

        nsnapshots = len(trajectory)
        H = trajectory[0].total_energy
        for snapshot_index in range(1, nsnapshots-1):
            H += trajectory[snapshot_index].kinetic_energy

        return H

    def computeActivity(self, trajectory):
        """
        Compute the activity of a given trajectory, defined in Ref. [1] as

        K[x(t)] = delta_t \sum_{t=0}^{t_obs} \sum_{j=1}^N [r_j(t+delta_t) - r_j(t)]^2

        RETURNS

        K (simtk.unit) - activity K[x(t)] for the specified trajectory

        """

        # Determine number of frames in trajectory.
        nframes = len(trajectory)

        # Compute activity of component A.
        K = 0.0 * self.delta_t * units.nanometers**2
        for frame_index in range(nframes-1):
            # Compute displacement of all atoms.
            delta_r = trajectory[frame_index+1].coordinates - trajectory[frame_index].coordinates
            # Compute contribution to activity K.
            K += self.delta_t * ((delta_r[0:self.N,:] / units.nanometers)**2).sum() * (units.nanometers**2)

        return K 

    def sampleTrajectory(self, trajectory):
        """
        Conduct one step of transition path sampling MCMC to generate a (correlated) trajectory sample.

        ARGUMENTS

        trajectory (Trajectory) - a previous trajectory sample from the TPS ensemble

        RETURN VALUES

        trajectory (Trajectory) - new sampled trajectory (correlated with previous trajectory sample)

        NOTE

        The new trajectory may share Snapshot objects from the old trajectory; modification of these objects
        will result in both old and new trajectories being updated.  Make a deep copy if necessary to keep these
        objects fully independent.

        """

        # Determine length of trajectory
        nframes = len(trajectory)

        # Compute value of activity K[x(t)]
        K_old = self.computeActivity(trajectory)

        # Choose a shooting or shift move
        SHOOT_PROBABILITY = 0.5 # probability of picking a shooting move or a shift move
        if (numpy.random.rand() < SHOOT_PROBABILITY):
            # Shoot part of a new trajectory
            # Pick a timeslice to shoot from.
            frame_index = numpy.random.random_integers(1, nframes-2)
            # Pick a shooting direction.
            if (numpy.random.rand() < 0.5):
                # Shoot forward.
                if self.verbose: print "Shooting forward from frame %d" % frame_index                
                partial_trajectory = self.generateTrajectory(trajectory[frame_index].coordinates, nframes - frame_index - 1)
                trial_trajectory = trajectory[0:frame_index] + partial_trajectory
            else:
                # Shoot backwards
                if self.verbose: print "Shooting backward from frame %d" % frame_index                                
                partial_trajectory = self.generateTrajectory(trajectory[frame_index].coordinates, frame_index)
                partial_trajectory.reverse()
                trial_trajectory = partial_trajectory[:-1] + trajectory[frame_index:]
        else:
            # Shift trajectory.
            # Pick a timeslice to form start of new trajectory.
            # Don't allow shifting by zero -- this screws with python indexing.
            nshift = numpy.random.random_integers(1, nframes-2)
            # Pick a shooting direction.
            if (numpy.random.rand() < 0.5):
                if self.verbose: print "Shifting by +%d" % nshift
                # Shoot forward from end.
                partial_trajectory = self.generateTrajectory(trajectory[-1].coordinates, nshift)
                trial_trajectory = trajectory[nshift:-1] + partial_trajectory
            else:
                # Shoot backwards from beginning.
                if self.verbose: print "Shifting by -%d" % nshift
                partial_trajectory = self.generateTrajectory(trajectory[0].coordinates, nshift)
                partial_trajectory.reverse()
                trial_trajectory = partial_trajectory[:-1] + trajectory[0:-nshift]

        # Compute new activity
        K_trial = self.computeActivity(trial_trajectory)

        # TODO: Test for crystallization
#        ncrystalframes = 0
#        for snapshot in trial_trajectory:
#            if is_crystal(snapshot.coordinates):
#                ncrystalframes += 1
#        reject = False
#        if ncrystalframes > threshhold:
#            reject = True
        
        # Accept or reject according to s-field.
        #print "s * (sigma**2 * delta_t)       = %.5f" % (self.s / self.s_reduced_unit)
        #print "K_old / (N * t_obs * sigma**2)   = %.5f" % (K_old / self.K_reduced_unit)
        #print "K_trial / (N * t_obs * sigma**2) = %.5f" % (K_trial / self.K_reduced_unit)
        log_P_accept = - self.s * (K_trial - K_old)
        #print "log_P_accept = %f" % log_P_accept
        self.nattempted += 1
        if (log_P_accept > 0.0) or (numpy.random.rand() < math.exp(log_P_accept)):
            # Accept trajectory move.
            if self.verbose: print "Accepted."
            self.naccepted += 1
            K_old = K_trial
            trajectory = trial_trajectory
        else:
            # Move was rejected
            if self.verbose: print "Rejected."            
            pass

        return trajectory

#=============================================================================================
# REPLICA-EXCHANGE TRANSITION PATH SAMPLING
#=============================================================================================

class ReplicaExchangeTPS(object):
    """
    Replica-exchange transition path sampling among different values of s-field.

    """

    def __init__(self, system, ensembles, trajectories, ncfilename, comm=None):
        """
        Initialize replica-exchange transition path sampling.

        ARGUMENTS

        system (simtk.chem.openmm.State) - system object
        ensemble (list of TransitionPathSampling objects)
        trajectories - trajectory or list of trajectories to initialize TPS with
        ncfilename - name of NetCDF file to create or resume from

        OPTIONAL ARGUMENTS

        comm - MPI communicator to run in parallel (default: None)
        

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
        self.comm = comm
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
            raise Error

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

        # Turn off verbosity if not master node.
        if self.comm is not None:
            if self.comm.rank != 0: self.verbose = False
  
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

        TODO

        * Create a parallel implementation when OpenMM supports this.

        """        
        
        # Propagate all replicas.
        if self.comm is None:
            # Serial version
            for replica_index in range(self.nstates):
                state_index = self.replica_states[replica_index]
                trajectory = self.trajectories[replica_index]
                trajectory = self.ensembles[state_index].sampleTrajectory(trajectory)
                self.trajectories[replica_index] = trajectory
        else:
            # Parallel version.
            print "Parallel _propagate_replicas: rank %d / %d" % (self.comm.rank, self.comm.size)

            # Share state information with all replicas.
            if self.comm.rank == 0: print "Sharing state information..."
            initial_time = time.time()            
            if self.comm.rank != 0: self.replica_states = None            
            self.replica_states = self.comm.bcast(self.replica_states, root=0)
            final_time = time.time()
            #self.comm.barrier()            
            if self.comm.rank == 0: print "Sharing state information: elapsed time %.3f s" % (final_time - initial_time)            

            # Run just this node's share of replicas.
            if self.comm.rank == 0: print "Running trajectories...."
            initial_time = time.time()
            for replica_index in range(self.comm.rank, self.nstates, self.comm.size):
                print "node %d / %d : running replica %d / %d" % (self.comm.rank, self.comm.size, replica_index, self.nstates)
                state_index = self.replica_states[replica_index]                
                trajectory = self.trajectories[replica_index]                
                trajectory = self.ensembles[state_index].sampleTrajectory(trajectory)
                self.trajectories[replica_index] = trajectory
            final_time = time.time()
            self.comm.barrier()
            if self.comm.rank == 0: print "Running trajectories: elapsed time %.3f s" % (final_time - initial_time)
            
            # Send trajectories back to root node.
            if self.comm.rank == 0: print "Synchronizing trajectories..."
            initial_time = time.time()
            gather = self.comm.allgather(self.trajectories[self.comm.rank:self.nstates:self.comm.size])
            for replica_index in range(self.nstates):
                source = replica_index % self.comm.size # node with trajectory data
                index = replica_index // self.comm.size # index within trajectory batch
                self.trajectories[replica_index] = gather[source][index]                    
            final_time = time.time()
            #self.comm.barrier()
            if self.comm.rank == 0: print "Synchronizing trajectories: elapsed time %.3f s" % (final_time - initial_time)

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
        for replica_index in range(self.nstates):
            K = self.activities[replica_index]
            H = self.path_hamiltonians[replica_index]
            for state_index in range(self.nstates):
                s = self.ensembles[state_index].s
                beta = self.ensembles[state_index].beta
                self.log_P_kl[replica_index,state_index] = - beta * H - s * K
        
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

        # Determine number of swaps to attempt to ensure thorough mixing.
        # TODO: Replace this with analytical result computed to guarantee sufficient mixing.
        nswap_attempts = self.nstates**4 # number of swaps to attempt
        
        # Allocate storage to keep track of mixing.
        Nij_proposed = numpy.zeros([self.nstates,self.nstates], numpy.float32) # Nij_proposed[i][j] is the number of swaps proposed between states i and j, prior of 1
        Nij_accepted = numpy.zeros([self.nstates,self.nstates], numpy.float32) # Nij_proposed[i][j] is the number of swaps proposed between states i and j, prior of 1

        # Show log P
        if self.verbose:
            print "log_P[replica,state] ="
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
        #ncvar_trajectory_velocities  = ncfile.createVariable('trajectory_velocities',  'f', ('replica','frame','atom','spatial'))
        ncvar_trajectory_potential   = ncfile.createVariable('trajectory_potential',   'f', ('replica','frame'))
        ncvar_trajectory_kinetic     = ncfile.createVariable('trajectory_kinetic',     'f', ('replica','frame'))

        ncvar_states      = ncfile.createVariable('states', 'i', ('iteration','replica'))
        ncvar_activities  = ncfile.createVariable('activities', 'f', ('iteration','replica'))
        ncvar_path_hamiltonians = ncfile.createVariable('path_hamiltonians', 'f', ('iteration','replica'))        
        ncvar_log_probabilities = ncfile.createVariable('log_probabilities', 'f', ('iteration','replica','replica'))
        ncvar_mixing      = ncfile.createVariable('mixing', 'f', ('iteration','replica','replica'))
        
        # Define units for variables.
        setattr(ncvar_trajectory_coordinates, 'units', 'nm')
        #setattr(ncvar_trajectory_velocities,  'units', 'nm/ps')
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

        if self.comm is not None:
            # Only master node writes in parallel runs.
            if self.comm.rank != 0: return        

        # DEBUG
        ensemble = self.ensembles[0]
        setattr(self.ncfile, 'sKfactor', ensemble.N * ensemble.t_obs / ensemble.delta_t) # the quantity by which s * K must be multipled to obtain log weight
        setattr(self.ncfile, 'betaHfactor', 1.0) # the quantity by which beta * H must be multipled to obtain log weight

        # Store trajectories.
        for replica_index in range(self.nstates):
            trajectory = self.trajectories[replica_index]
            for frame_index in range(self.nframes):                
                self.ncfile.variables['trajectory_coordinates'][replica_index,frame_index,:,:] = (trajectory[frame_index].coordinates / units.nanometers).astype(numpy.float32)
                #self.ncfile.variables['trajectory_velocities'][replica_index,frame_index,:,:] = (trajectory[frame_index].velocities / (units.nanometers / units.picoseconds)).astype(numpy.float32)
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
        self.ncfile.variables['log_probabilities'][self.iteration,:,:] = self.log_P_kl[:,:]
        
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
                #v = ncfile.variables['trajectory_velocities'][replica_index,frame_index,:,:].astype(numpy.float32).copy()
                #velocities = units.Quantity(v, units.nanometers / units.picoseconds)                
                V = ncfile.variables['trajectory_potential'][replica_index,frame_index]
                potential_energy = units.Quantity(V, units.kilojoules_per_mole)
                T = ncfile.variables['trajectory_kinetic'][replica_index,frame_index]
                kinetic_energy = units.Quantity(T, units.kilojoules_per_mole)
                #snapshot = Snapshot(coordinates=coordinates, velocities=velocities, kinetic_energy=kinetic_energy, potential_energy=potential_energy)
                snapshot = Snapshot(coordinates=coordinates, kinetic_energy=kinetic_energy, potential_energy=potential_energy)                
                trajectory.append(snapshot)
            self.trajectories.append(trajectory)

        # Restore state information.
        self.replica_states = ncfile.variables['states'][self.iteration,:].copy()

        # Restore log probabilities.
        self.log_P_kl = ncfile.variables['log_probabilities'][self.iteration,:,:]
        
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
        if (self.comm is None) or (self.comm.rank == 0):
            self.ncfile = netcdf.Dataset(self.store_filename, 'a') # for netCDF4
        else:
            self.ncfile = None
        
        return

#=============================================================================================
# MAIN AND TESTS
#=============================================================================================

class Driver(object):
    def __init__(self):
        pass

    def run(self):
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
        platform = openmm.Platform.getPlatformByName('OpenCL') # fast GPU platform
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
#        # Create a Context.
#        timestep = 1.0 * units.femtoseconds
#        integrator = openmm.VerletIntegrator(timestep)
#        context = openmm.Context(system, integrator, platform)
#        # Set coordinates.
#        context.setPositions(coordinates)
#        # Compute initial energy.
#        state = context.getState(getEnergy=True)
#        initial_potential = state.getPotentialEnergy()
#        # Minimize.
#        openmm.LocalEnergyMinimizer.minimize(context)    
#        # Compute final energy.
#        state = context.getState(getEnergy=True, getCoordinates=True)
#        final_potential = state.getPotentialEnergy()
#        coordinates = state.getPositions(asNumpy=True)        
            # Minimize the initial coordinates.
            minimizer = optimize.LBFGSMinimizer(system, verbose=True, platform=platform)
            coordinates = minimizer.minimize(coordinates)
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
        ensembles = [ TransitionPathSampling(system, N, NA, mass, epsilon, sigma, timestep, nsteps_per_frame, nframes, quenched_temperature, s * s_reduced_unit, platform=platform) for s in svalues ]

        trajectory = None
        if seed:
            # Generate an initial trajectory at zero field
            print "Generating seed trajectory for TPS..."
            trajectory = ensembles[0].generateTrajectory(coordinates, nframes)

        # Initialize replica-exchange TPS simulation.
        print "Initializing replica-exchange TPS..."
        simulation = ReplicaExchangeTPS(system, ensembles, trajectory, ncfilename)

        # Set simulation parameters.
        simulation.number_of_iterations = 10000

        # Run simulation
        print "Running replica-exchange TPS..."
        simulation.run()

        print "Work Done :) Go home and have a beer."
    
#=============================================================================================
# This driver sets up a simulation for multiple temperatures and s-values.
#=============================================================================================
class MultitempDriver(object):
    def __init__(self):
        pass

    def run(self):
        print "multitemp driver"
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

        # Set up device to bind to.
        print "MPI..."
        from mpi4py import MPI # MPI wrapper    
        hostname = os.uname()[1]
        ngpus = 2
        comm = MPI.COMM_WORLD # MPI communicator
        deviceid = comm.rank % ngpus
        platform = openmm.Platform.getPlatformByName('OpenCL') # fast GPU platform
        platform.setPropertyDefaultValue('OpenCLDeviceIndex', '%d' % deviceid) # use first OpenCL device
        print "node '%s' deviceid %d / %d, MPI rank %d / %d" % (hostname, deviceid, ngpus, comm.rank, comm.size)
        # Make sure random number generators have different seeds.
        seed = numpy.random.randint(sys.maxint - comm.size) + comm.rank
        numpy.random.seed(seed)

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

        # Create a number of transition path sampling ensembles at different values of the field parameter s and different temperatures
        s_reduced_unit = 1.0 / (sigma**2 * delta_t) # reduced units for s field
        T_reduced_unit = epsilon / kB # reduced units for temperature

        svalues = [0.00, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, -0.005]
        Tvalues = [quenched_temperature/T_reduced_unit for s in svalues] # use quenched temperature for all states

        print "Initializing transition path sampling ensembles..."
        ensembles = [ TransitionPathSampling(system, N, NA, mass, epsilon, sigma, timestep, nsteps_per_frame, nframes, T * T_reduced_unit, s * s_reduced_unit, platform=platform) for (s,T) in zip(svalues,Tvalues) ]

        trajectory = None
        if seed:
            # Generate an initial trajectory at zero field
            print "Generating seed trajectory for TPS..."
            trajectory = ensembles[0].generateTrajectory(coordinates, nframes)

        # Initialize replica-exchange TPS simulation.
        print "Initializing replica-exchange TPS..."
        simulation = ReplicaExchangeTPS(system, ensembles, trajectory, ncfilename, comm=comm)

        # Set simulation parameters.
        simulation.number_of_iterations = 20000

        # Run simulation
        print "Running replica-exchange TPS..."
        print "(I'll buy you a beer if you actually get this to work.)"
        simulation.run()

        print "And there was much rejoicing!"

#if __name__ == "__main__":
# Run new-and-improved multiple temperature and s value driver.
driver = MultitempDriver()
driver.run()
    

    
   
  
