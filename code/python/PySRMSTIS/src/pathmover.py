'''
Created on 19.07.2014

@author: jan-hendrikprinz
'''

import math
from trajectory import Trajectory
import logger


class PathMover(object):
    def __init__(self):
        
        # A list of stopping criterions
        self.stopper = None
        # A a ShootingPointSelector
        self.selector = None
        pass
    
    def move(self, trajectory):
        shooting_point, bias = self.selector.select(path)
        self._generate()
        self._accept()
        
        return old, new, accept
        pass

class ForwardShoot(PathMover):

    def move(self, trajectory):
        # Determine length of trajectory
        nframes = len(trajectory)
        xframes = self.simulator.n_frames_max
        
        t = Trajectory()
        
        # Compute ptah Hamiltonian
        l_old = 1.0 * len(trajectory)
        max_length_prop = math.ceil(l_old / np.random.rand())

        if max_length_prop < xframes:
            xframes = max_length_prop

        # Stopper
        stopper = self.stopper
        selector = self.selector
        acceptor = self.acceptor
        
        # Choose a shooting or shift move
        # Shoot part of a new trajectory
        # Pick a timeslice to shoot from.
        # TODO: This could be changed to more transition regions
        frame_index = selector.pick(trajectory)
        # Pick a shooting direction.
        print "Shooting forward from frame %d" % frame_index
        l_max = xframes - frame_index - 1

        partial_trajectory = self.simulator.generate(trajectory[start_frame], l_max, stopping = stopper)
        print "Trial was ", len(partial_trajectory) + frame_index, " long"                

        if len(partial_trajectory) == l_max + 1:
            print "Rejected. Too long"            
            trial_trajectory = trajectory
        else:
            print "Accepted."
            trial_trajectory = trajectory[0:frame_index] + partial_trajectory

        return trial_trajectory

class BackwardShoot(PathMover):
    def move(self, trajectory):
        
                # Shoot backwards
                print "Shooting backward from frame %d" % frame_index     
                l_max = xframes - nframes + frame_index
                partial_trajectory = self.simulator.generate(trajectory[frame_index], l_max, stopping = stopper)
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
                partial_trajectory = self.simulator.generate(trajectory[-1], xframes - nframes + nshift, stopping = stopper)
                trial_trajectory = trajectory[nshift:-1] + partial_trajectory
            else:
                # Shoot backwards from beginning.
                print "Shifting by -%d" % nshift
                partial_trajectory = self.simulator.generate(trajectory[0], xframes - nframes + nshift, stopping = stopper)
                partial_trajectory.reverse()
                trial_trajectory = partial_trajectory[:-1] + trajectory[0:-nshift]

        # Compute new path Hamiltonian
#        H_trial = trial_trajectory.pathHamiltonian()

        
#        log_P_accept = - self.beta * (H_trial - H_old)
        #print "log_P_accept = %f" % log_P_accept
        


class ReplicaExchange(object):
    def do_move(self, allpaths, state):
        pass

class MinusMove(object):
    def do_move(self, allpaths, state):
        pass

class PathReversal(object):
    def do_move(self, allpaths, state):
        pass
