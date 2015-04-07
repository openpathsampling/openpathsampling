'''
@author: JD Chodera
@author: ASJS Mey
@author: JH Prinz
'''

import numpy as np 
import math
from simtk.unit import nanosecond, picosecond, nanometers, nanometer, picoseconds, femtoseconds, femtosecond

import simtk.unit as units


from trajectory import Trajectory

kB = units.BOLTZMANN_CONSTANT_kB * units.AVOGADRO_CONSTANT_NA

class TransitionInterfaceSampling(object):
    """
    General Transition path sampling with stopping criterion

    """

    def __init__(self, simulator, interfaces):
        """
        Initialize a transition path sampling simulation for a given simulation context

        ARGUMENTS

        context (simtk.openmm) - the openmm simulation context

        OPTIONAL ARGUMENTS

        platform (simtk.chem.openPlatform) - platform to use for OpenMM pathsimulators
        
        """

        # Store reference to the openmm simulation system.
        self.simulator = simulator
        self.interfaces = interfaces
        self.beta = 1/ (298 * units.kelvin * kB)
        
        print "beta : ", self.beta

        return

    def sampleTrajectory(self, trajectory):
        """
        Conduct one step of transition path sampling MCMC to generate a (correlated) trajectory sample.

        ARGUMENTS

        trajectory (Trajectory) - a previous trajectory sample from the TPS ensemble

        RETURN VALUES

        trajectory (Trajectory) - new sampled trajectory (correlated with previous trajectory sample)

        NOTE

        Snapshot objects should be immutable!!!! At least where we can make sure of this.
        
        Does it make sense to roll the probability to reject a trajectory before hand and set the maximal length accordingly.
        If we will reject a too long trajectory anyway we can stop once it became too long, right?
        
        """

        # Determine length of trajectory
        nframes = len(trajectory)
        xframes = self.simulator.n_frames_max
        
        t = Trajectory()
        t.pathHamiltonian()
        
        # Compute ptah Hamiltonian
#        H_old = trajectory.pathHamiltonian()
        l_old = 1.0 * len(trajectory)
        max_length_prop = math.ceil(l_old / np.random.rand())
        if max_length_prop < xframes:
            xframes = max_length_prop
            
        # Stopper
        stopper = self.interfaces.in_core
        
        print "Using trajectory of length : ", nframes, " (max allowed : ", xframes, " )",
        
        # Choose a shooting or shift move
        SHOOT_PROBABILITY = 1.0 # probability of picking a shooting move or a shift move
        if (np.random.rand() < SHOOT_PROBABILITY):
            # Shoot part of a new trajectory
            # Pick a timeslice to shoot from.
            # TODO: This could be changed to more transition regions
            frame_index = np.random.random_integers(1, nframes-2)
            # Pick a shooting direction.
            if (np.random.rand() < 0.5):
                # Shoot forward.
                print "Shooting forward from frame %d" % frame_index
                l_max = xframes - frame_index - 1
                partial_trajectory = self.simulator.generate(trajectory[frame_index], l_max, running = stopper)
                print "Trial was ", len(partial_trajectory) + frame_index, " long"                
                if len(partial_trajectory) == l_max + 1:
                    print "Rejected. Too long"            
                    trial_trajectory = trajectory
                else:
                    print "Accepted."
                    trial_trajectory = trajectory[0:frame_index] + partial_trajectory
            else:
                # Shoot backwards
                print "Shooting backward from frame %d" % frame_index     
                l_max = xframes - nframes + frame_index
                partial_trajectory = self.simulator.generate(trajectory[frame_index], l_max, running = stopper)
                print "Trial was ", len(partial_trajectory) + nframes - frame_index, " long"                
                
                if len(partial_trajectory) == l_max + 1:
                    print "Rejected. Too long"            
                    trial_trajectory = trajectory
                else:
                    print "Accepted."
                    partial_trajectory.reverse()
                    trial_trajectory = partial_trajectory[:-1] + trajectory[frame_index:]
        else:
            # Shift trajectory.
            # Pick a timeslice to form start of new trajectory.
            # Don't allow shifting by zero -- this screws with python indexing.
            nshift = np.random.random_integers(1, nframes-2)
            # Pick a shooting direction.
            if (np.random.rand() < 0.5):
                print "Shifting by +%d" % nshift
                # Shoot forward from end.
                partial_trajectory = self.simulator.generate(trajectory[-1], xframes - nframes + nshift, running = stopper)
                trial_trajectory = trajectory[nshift:-1] + partial_trajectory
            else:
                # Shoot backwards from beginning.
                print "Shifting by -%d" % nshift
                partial_trajectory = self.simulator.generate(trajectory[0], xframes - nframes + nshift, running = stopper)
                partial_trajectory.reverse()
                trial_trajectory = partial_trajectory[:-1] + trajectory[0:-nshift]

        # Compute new path Hamiltonian
#        H_trial = trial_trajectory.pathHamiltonian()

        
#        log_P_accept = - self.beta * (H_trial - H_old)
        #print "log_P_accept = %f" % log_P_accept
        

        return trial_trajectory
    


class TransitionInterfaceSampling_David(object):

    def __init__(self):
        movers = []

    def run(self):
        for step in range(nsteps):
            movers[rand_num].do_move(allpaths, state)

    
# Store interface definitions -> might go to a class like collective variable or lambda definition
# Store attempts with initial trajectory, final trajectory, result accepted/rejected wrong ensemble, rejected length, etc..., also which ensemble is sampled from
# Store a list of which trajectory belongs to which ensemble. Hope this is enough to also add all rejections to an ensemble later on.