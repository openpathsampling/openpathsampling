'''
Created on 19.07.2014

@author: Jan-Hendrik Prinz
'''

from shooting import ShootingPoint
from ensemble import ForwardAppendedTrajectoryEnsemble, BackwardPrependedTrajectoryEnsemble, LengthEnsemble

import numpy as np

class PathMover(object):
    '''
    A PathMover is the description of how to generate a new path from an old one.
    
    Notes
    -----
    
    Basically this describes the proposal step for a MC in path space.
    
    We might detach this from the acceptance step?!?!?
    This would mean that a PathMover needs only an old trajectory and gives a new one.
    
    For example a ForwardShoot then uses a shooting point selector and runs a new trajectory and combine them to get
    a new one.
    
    After the move has been made, we can retrieve information about the move, as well as the new trajectory from the
    PathMover object
    
    Attributes
    ----------
    simulator : Simulator
        the attached simulator used to generate new trajectories
        
    '''
    
    simulator = None        
    
    def __init__(self):
        
        # An ensemble that is at the same time triggering the stopping criterion and the final acceptance and the end.
        # This is because the goal is to sample trajectories in a specific ensemble. So we want to generate and stop
        # as soon as this cannot be fulfilled anymore. Some of the conditions cannot be checked during runtime, so we
        # have to do that at the end to make sure.
        
        self.ensemble = None
        
        self.start = None
        self.final = None        
        pass
        
    def move(self, trajectory):
        '''
        Run the generation starting with the initial trajectory specified.

        Parameters
        ----------
        trajectory : Trajectory
            the initial trajectory
        
        Returns
        -------        
        final : Trajectory()
            the final new proposed trajectory
        
        Notes
        -----
        After this command additional information can be accessed from this object
        '''
        self.start = trajectory
        self.final = trajectory
                        
        return self.final

    @property
    def proposal_bias(self):
        '''
        Return the proposal bias necessary to correct for an asymmetric proposal.
        
        Notes
        -----
        This is effectively the ratio of proposal probabilities for a mover. For symmetric proposal
        this is one. In the case of e.g. Shooters this depends on the used ShootingPointSelector and
        the start and final trajectory.
        
        I am not sure if it makes sense that to define it this way, but for Shooters this is,
        what we need for the acceptance step in addition to the check if we have a trajectory of
        the target ensemble.

        What about Minus Move and PathReversalMove?
        '''
        return 1.0

class ShootMover(PathMover):
    '''
    A pathmover that implements a general shooting algorithm that generates a sample from a specified ensemble
    '''    
    def __init__(self, selector, ensemble):
        super(ShootMover, self).__init__()
        self.selector = selector
        self.ensemble = ensemble
        self.length_stopper = LengthEnsemble(0)
    
    @property
    def generator_stopper(self):
        '''
        the ensemble used by the generator 
        
        Notes
        -----
        contains the status of the generator if there has been a too long trajectory generated
        
        '''
        return PathMover.simulator.max_length_stopper
    
    @property
    def proposal_bias(self):
        '''
        Return the proposal bias for Shooting Moves. These are given by the ratio of partition functions
        '''
        return self.start_point.Z / self.final_point.Z
    
    def _generate(self):
        self.final = self.start
    
    def move(self, trajectory):        
        self.start = trajectory        
        self.start_point = self.selector.pick(self.start)
                        
        max_frames = self.simulator.n_frames_max - trajectory.frames        
        self.length_stopper.length = max_frames
        
        self._generate()

        self.accepted = self.ensemble(self.final)
        
        if self.accepted:
            rand = np.random.random()
            print 'Proposal probability', self.proposal_bias, '/ random :', rand
            if (rand < self.proposal_bias):
                return self.final
            
        return self.start
    
    
class ForwardShootMover(ShootMover):    
    '''
    A pathmover that implements the forward shooting algorithm
    '''
    def _generate(self):
        print "Shooting forward from frame %d" % self.start_point.index
        
        # Run until one of the stoppers is triggered
        partial_trajectory = PathMover.simulator.generate(
                                     self.start_point.snapshot, 
                                     running = [ForwardAppendedTrajectoryEnsemble(
                                                      self.ensemble, 
                                                      self.start[0:self.start_point.index]
                                                ).forward, self.length_stopper.forward]
                                     )        
        self.final = self.start[0:self.start_point.index] + partial_trajectory    
        self.final_point = ShootingPoint(self.selector, self.final, self.start_point.index)

        pass
    
class BackwardShootMover(ShootMover):    
    '''
    A pathmover that implements the backward shooting algorithm
    '''
    def _generate(self):
        print "Shooting backward from frame %d" % self.start_point.index
#        print self.start_point
        # Run until one of the stoppers is triggered
        partial_trajectory = PathMover.simulator.generate(
                                     self.start_point.snapshot.reversed_copy(), 
                                     running = [BackwardPrependedTrajectoryEnsemble(
                                                     self.ensemble, 
                                                     self.start[self.start_point.index + 1:]
                                                ).backward, self.length_stopper.backward]
                                     )
        
        
        self.final = partial_trajectory.reversed + self.start[self.start_point.index + 1:]    
        self.final_point = ShootingPoint(self.selector, self.final, partial_trajectory.frames - 1)

#        print self.start_point.snapshot.velocities
#        print self.start_point.snapshot.momentum.idx
#        print self.final_point.snapshot.velocities
#        print self.final_point.snapshot.momentum.idx

        
        pass

class MixedMover(PathMover):
    '''
    Defines a mover that picks a over from a set of movers with specific weights.
    
    Notes
    -----
    Channel functions from self.mover to self. Does NOT work yet. Think about a good way to implement this...
    '''
    def __init__(self, movers, weights = None):
        self.movers = movers        
        self.mover = None

        if weights is None:
            self.weights = [1.0] * len(movers)
        else:
            self.weights = weights
        pass
    
    def move(self, trajectory):
        rand = np.random.random() * sum(self.weights)
        idx = 0
        prob = self.weights[0]
        while prob <= rand and idx < len(self.weights):
            idx += 1
            prob += self.weights[idx]
            
        self.mover = self.movers[idx]

        return self.mover.move(trajectory)
    
    def __getattr__(self, attr):
        # relay the attributes from the selected mover to the class itself for easier access
        return getattr(self.mover, attr)


#############################################################
# The following moves still need to be implemented. Check what excactly they do
#############################################################

class MinusMove(object):
    def do_move(self, allpaths, state):
        pass

class PathReversal(object):
    def do_move(self, allpaths, state):
        pass


#############################################################
# The following move should be moved to RETIS and just uses moves. It is not a move itself
#############################################################

class ReplicaExchange(object):
    def do_move(self, allpaths, state):
        pass