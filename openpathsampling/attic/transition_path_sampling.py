'''
@author: JD Chodera
@author: ASJS Mey
@author: JH Prinz
'''

from simtk.unit import nanosecond, picosecond, nanometers, nanometer, picoseconds, femtoseconds, femtosecond


from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *

from trajectory import Trajectory
from snapshot import Snapshot

from openmmtools.integrators import VVVRIntegrator


#=============================================================================================
# GLOBAL CONSTANTS
#=============================================================================================

kB = BOLTZMANN_CONSTANT_kB * AVOGADRO_CONSTANT_NA

#=============================================================================================
# TRANSITION PATH SAMPLING
#=============================================================================================

class TransitionPathSampling(object):
    """
    General Transition path sampling with stopping criterion

    """

    def __init__(self, system, platform=None):
        """
        Initialize a transition path sampling simulation with activity field strength s.

        ARGUMENTS

        system (simtk.chem.openSystem) - the molecular mechanics system
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

        platform (simtk.chem.openPlatform) - platform to use for OpenMM calculations
        
        """

        # Store local copy of System.
        self.system = system

#         self.mass = mass
#         self.epsilon = epsilon
#         self.sigma = sigma
# 
#         self.timestep = timestep
#         self.nsteps_per_frame = nsteps_per_frame
#         self.delta_t = timestep * nsteps_per_frame
# 

        # Compute thermal energy and inverse temperature from specified temperature.
        self.kT = kB * self.temperature # thermal energy
        self.beta = 1.0 / self.kT # inverse temperature

        # Create a Context for integration.
        if platform:
            self.context = Context(self.system, self.integrator, platform)
        else:
            self.context = Context(self.system, self.integrator)            

        # Store reduced units
#        self.t_obs = nframes * self.delta_t
#        self.s_reduced_unit = 1.0 /  (self.sigma**2 * self.delta_t)
#        self.K_reduced_unit = (self.N * self.t_obs * self.sigma**2)
#        self.H_reduced_unit = self.epsilon # reduced unit for path Hamiltonian (energy)
#        self.beta_reduced_unit = 1.0 / self.epsilon # reduced unit for inverse temperature

        self.temperature = 298.0 * kelvin
        self.collision_rate = 91.0 / picoseconds
        self.timestep = 1.0 * femtoseconds
        self.integrator = VVVRIntegrator(self.temperature, self.collision_rate, self.timestep)

        return

    def generateTrajectory(self, x0, nframes):
        """
        Generate a velocity Verlet trajectory consisting of ntau segments of tau_steps in between storage of Snapshots and randomization of velocities.

        ARGUMENTS
        
        x0 (coordinate set) - initial coordinates; velocities will be assigned from Maxwell-Boltzmann distribution
        nframes (int) - number of trajectory segments to generate

        RETURNS

        trajectory (list of Snapshot) - generated trajectory of initial conditions, including initial coordinate set

        NOTES
        
        This exists in OpenMM. Check with John on how to best get snapshots
        This routine generates a velocity Verlet trajectory for systems without constraints by wrapping the OpenMM 'VerletIntegrator' in two half-kicks of the velocity.
        
        """

        # Set initial positions
        self.context.setPositions(x0)

        # Store initial state for each trajectory segment in trajectory.
        trajectory = Trajectory()

        # Construct mass vector.
        nparticles = self.system.getNumParticles()
        mass = Quantity(numpy.zeros([nparticles,3], numpy.float64), amu)
        for particle_index in range(nparticles):
            mass[particle_index,:] = self.system.getParticleMass(particle_index)

            
        # Assign velocities from Maxwell-Boltzmann distribution            
        self.context.setVelocitiesToTemperature(self.temperature)
        
        # Store initial snapshot of trajectory segment.
        snapshot = Snapshot(context=self.context)
        trajectory.forward(snapshot)
        
        # Propagate dynamics by velocity Verlet.
        in_void = True
        
        self.frame = 0
        
        while self.frame < self.xframes and in_void == True:
            
            # Do integrator x steps
            self.integrator.step()
            
            
            
            # Check if reached a core set
            in_void = self.check_void()
            

        # Store final snapshot of trajectory.
        snapshot = Snapshot(self.context)
        trajectory.forward(snapshot)

        return trajectory

    def logEquilibriumTrajectoryProbability(self, trajectory):
        """
        Compute the log equilibrium probability (up to an unknown additive constant) of an unbiased trajectory evolved according to Verlet dynamics with Andersen thermostatting.

        ARGUMENTS

        trajectory (Trajectory) - the trajectory

        RETURNS

        log_q (float) - the log equilibrium probability of the trajectory
        
        NOTES
        This might be better places into the trajectory class. The trajectory should know the system and ensemble? and so it is not necessarily 
        TPS specific

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
        
        NOTES
        Could also more efficiently placed in the trajectory class

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
        
        NOTES
        
        Place also into trajectory class

        """

        # Determine number of frames in trajectory.
        nframes = len(trajectory)

        # Compute activity of component A.
        K = 0.0 * self.delta_t * nanometers**2
        for frame_index in range(nframes-1):
            # Compute displacement of all atoms.
            delta_r = trajectory[frame_index+1].coordinates - trajectory[frame_index].coordinates
            # Compute contribution to activity K.
            K += self.delta_t * ((delta_r[0:self.N,:] / nanometers)**2).sum() * (nanometers**2)

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
        
        Snapshot objects should be immutable!!!! At least where we can make sure of this.
        
        This should be changed to TIS, but should maybe keep all the trajetory idea, etc and just extend to TIS

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
            # TODO: This could be changed to more transition regions
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
        if (log_P_accept > 0.0) or (numpy.random.rand() < math.exp(log_P_accept)):
            # Accept trajectory move.
            print "Accepted."
            self.naccepted += 1
            K_old = K_trial
            trajectory = trial_trajectory
        else:
            # Move was rejected
            print "Rejected."            
            pass

        return trajectory