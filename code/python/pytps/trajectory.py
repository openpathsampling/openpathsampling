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

import UserList

#=============================================================================================
# SIMULATION SNAPSHOT 
#=============================================================================================

class Snapshot(object):
    """
    Simulation snapshot.

    """
    def __init__(self, context=None, coordinates=None, velocities=None, box_vectors=None, potential_energy=None, kinetic_energy=None):
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
            self.velocities = state.getVelocities(asNumpy=True)
            self.box_vectors = state.getPeriodicBoxVectors() # TODO: set asNumpy=True once bug in OpenMM is fixed
            self.potential_energy = state.getPotentialEnergy()
            self.kinetic_energy = state.getKineticEnergy()
        else:
            import copy
            if coordinates is not None: self.coordinates = copy.deepcopy(coordinates)
            if velocities is not None: self.velocities = copy.deepcopy(velocities)
            if box_vectors is not None: self.box_vectors = copy.deepcopy(box_vectors)
            if potential_energy is not None: self.potential_energy = copy.deepcopy(potential_energy)
            if kinetic_energy is not None: self.kinetic_energy = copy.deepcopy(kinetic_energy)                       

        # Check for nans in coordinates, and raise an exception if something is wrong.
        import numpy
        if numpy.any(numpy.isnan(self.coordinates)):
            raise Exception("Some coordinates became 'nan'; simulation is unstable or buggy.")

        return

    @property
    def total_energy(self):
        return self.kinetic_energy + self.potential_energy

#=============================================================================================
# SIMULATION TRAJECTORY
#=============================================================================================

class Trajectory(UserList.UserList):
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
        UserList.UserList.__init__(self)
        
        if trajectory is not None:
            # Try to make a copy out of whatever container we were provided
            for snapshot in trajectory:
                import copy
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
        UserList.UserList.reverse(self)

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

    def __init__(self, system, N, NA, mass, epsilon, sigma, timestep, nsteps_per_frame, nframes, temperature, s, platform_name=None):
        """
        Initialize a transition path sampling simulation with activity field strength s.

        ARGUMENTS

        system (pyopenmm.System) - the molecular mechanics system
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

        platform_name (string) - name of platform to run on, or None to select default
        
        """

        import numpy
        import simtk.unit as units

        # Store local copy of System.
        self.system = system

        self.mass = mass
        self.epsilon = epsilon
        self.sigma = sigma

        self.timestep = timestep
        self.nsteps_per_frame = nsteps_per_frame
        self.delta_t = timestep * nsteps_per_frame

        self.temperature = temperature

        # Compute thermal energy and inverse temperature from specified temperature.
        kB = units.BOLTZMANN_CONSTANT_kB * units.AVOGADRO_CONSTANT_NA        
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

        # Store reduced units
        self.t_obs = nframes * self.delta_t
        self.s_reduced_unit = 1.0 /  (self.sigma**2 * self.delta_t)
        self.K_reduced_unit = (self.N * self.t_obs * self.sigma**2)
        self.H_reduced_unit = self.epsilon # reduced unit for path Hamiltonian (energy)
        self.beta_reduced_unit = 1.0 / self.epsilon # reduced unit for inverse temperature

        # Store requested OpenMM platform name.
        self.platform_name = platform_name
        self.deviceid = 0

        # Form vectors of masses and sqrt(kT/m) for force propagation and velocity randomization.
        nparticles = system.getNumParticles()
        self.mass = units.Quantity(numpy.zeros([nparticles,3], numpy.float64), units.amu)
        for particle_index in range(nparticles):
            self.mass[particle_index,:] = self.system.getParticleMass(particle_index)
        kT = kB * temperature # thermal energy    
        self.sqrt_kT_over_m = units.Quantity(numpy.zeros([nparticles,3], numpy.float64), units.nanometers / units.picosecond)
        for particle_index in range(nparticles):
            self.sqrt_kT_over_m[particle_index,:] = units.sqrt(kT / self.mass[particle_index,0]) # standard deviation of velocity distribution for each coordinate for this atom

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
        import numpy
        velocities = self.sqrt_kT_over_m * numpy.random.standard_normal(size=(nparticles,3))

        if remove_com_velocity:
            import simtk.unit as units
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

        import pyopenmm

        # Create a Verlet integrator.
        integrator = pyopenmm.VerletIntegrator(self.timestep)

        # Create a Context for integration.
        if self.platform_name is not None:
            platform = pyopenmm.Platform.getPlatformByName(self.platform_name)
            platform.setPropertyDefaultValue('OpenCLDeviceIndex', '%d' % self.deviceid) # use first OpenCL device
            print "using OpenCL device %d" % self.deviceid
            context = pyopenmm.Context(self.system, integrator, platform)
        else:
            context = pyopenmm.Context(self.system, integrator)            

        # Set initial positions
        context.setPositions(x0)

        # Store initial state for each trajectory segment in trajectory.
        trajectory = Trajectory()

        # Generate trajectory segments.
        for frame_index in range(nframes):
            # Assign velocities from Maxwell-Boltzmann distribution
            velocities = self.assignMaxwellBoltzmannVelocities(remove_com_velocity=True)
            context.setVelocities(velocities)            

            # Store initial snapshot of trajectory segment.
            snapshot = Snapshot(context)
            trajectory.append(snapshot)
            
            # Propagate dynamics by velocity Verlet.
            # We only have leapfrog integrator available, so we wrap it in two half-kicks.
            # Back-kick by half a timestep to get ready for leapfrog integration.
            state = context.getState(getForces=True, getVelocities=True)
            force = state.getForces(asNumpy=True)
            velocities = state.getVelocities(asNumpy=True)
            velocities -= 0.5 * force/self.mass * self.timestep
            context.setVelocities(velocities)
            # Step using leapfrog.
            integrator.step(self.nsteps_per_frame)
            # Forward-kick by half a timestep to bring velocities into sync with positions.
            state = context.getState(getForces=True, getVelocities=True)
            force = state.getForces(asNumpy=True)
            velocities = state.getVelocities(asNumpy=True)
            velocities += 0.5 * force/self.mass * self.timestep            
            context.setVelocities(velocities)

        # Store final snapshot of trajectory.
        snapshot = Snapshot(context)
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
        import simtk.unit as units
        K = 0.0 * self.delta_t * units.nanometers**2
        for frame_index in range(nframes-1):
            # Compute displacement of all atoms.
            delta_r = trajectory[frame_index+1].coordinates - trajectory[frame_index].coordinates
            # Compute contribution to activity K.
            K += self.delta_t * ((delta_r[0:self.N,:] / units.nanometers)**2).sum() * (units.nanometers**2)

        return K 

    def sampleTrajectory(self, trajectory, deviceid=None):
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

        import time
        initial_time = time.time()
        if deviceid is not None:
            self.deviceid = deviceid
        
        # Determine length of trajectory
        nframes = len(trajectory)

        # Compute value of activity K[x(t)]
        K_old = self.computeActivity(trajectory)

        # Choose a shooting or shift move
        import numpy
        SHOOT_PROBABILITY = 0.5 # probability of picking a shooting move or a shift move
        if (numpy.random.rand() < SHOOT_PROBABILITY):
            # Shoot part of a new trajectory
            # Pick a timeslice to shoot from.
            frame_index = numpy.random.random_integers(1, nframes-2)
            # Pick a shooting direction.
            if (numpy.random.rand() < 0.5):
                # Shoot forward.
                print "Shooting forward from frame %d" % frame_index                
                partial_trajectory = self.generateTrajectory(trajectory[frame_index].coordinates, nframes - frame_index - 1)
                trial_trajectory = trajectory[0:frame_index] + partial_trajectory
            else:
                # Shoot backwards
                print "Shooting backward from frame %d" % frame_index                                
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
                print "Shifting by +%d" % nshift
                # Shoot forward from end.
                partial_trajectory = self.generateTrajectory(trajectory[-1].coordinates, nshift)
                trial_trajectory = trajectory[nshift:-1] + partial_trajectory
            else:
                # Shoot backwards from beginning.
                print "Shifting by -%d" % nshift
                partial_trajectory = self.generateTrajectory(trajectory[0].coordinates, nshift)
                partial_trajectory.reverse()
                trial_trajectory = partial_trajectory[:-1] + trajectory[0:-nshift]

        # Compute new activity
        K_trial = self.computeActivity(trial_trajectory)
        
        # Accept or reject according to s-field.        
        #print "s * (sigma**2 * delta_t)       = %.5f" % (self.s / self.s_reduced_unit)
        #print "K_old / (N * t_obs * sigma**2)   = %.5f" % (K_old / self.K_reduced_unit)
        #print "K_trial / (N * t_obs * sigma**2) = %.5f" % (K_trial / self.K_reduced_unit)
        log_P_accept = - self.s * (K_trial - K_old)
        #print "log_P_accept = %f" % log_P_accept
        self.nattempted += 1
        if (log_P_accept > 0.0) or (numpy.random.rand() < numpy.exp(log_P_accept)):
            # Accept trajectory move.
            print "Accepted."
            self.naccepted += 1
            K_old = K_trial
            trajectory = trial_trajectory
        else:
            # Move was rejected
            print "Rejected."            
            pass

        final_time = time.time()
        print "%.3f s elapsed" % (final_time - initial_time)

        return trajectory

    def test(self, trajectory):
        """
        Test function.
        """
        print "Testing."
        return 1

