'''
Created on 06.09.2014

@author: jan-hendrikprinz
'''

import numpy as np

class PathAlgorithm(object):
    '''
    PathAlgorithm contains a set of allowed PathMovers with weights, and an acceptance criterion
    '''


    def __init__(self, mover, acceptor):
        '''
        Constructor
        '''
        
        self.mover = mover
        self.acceptor = acceptor
        
        self.accepted = None
        
    def trial(self, trajectory):
        '''
        Run one of the PathMover object and then the acceptor to get a new path
        '''

        self.initial_traj = trajectory
        # Create a new trajectory                
        self.proposal_traj = self.mover.move(trajectory)

        # Compute biases from proposal
        self.proposal_bias = self.mover.selection_probability_ratio
        
        # and acceptor, which is usually zero or one
        self.acceptor_bias = self.acceptor.test(self.proposal_traj)
        
        self.prob = np.random.random()
        
        if self.prob < self.proposal_bias * self.acceptor_bias:
            # Accept
            self.accepted = True
            
        else:
            # Reject
            self.accepted = False
        
        return self.next()
    
    def next(self):
        if self.accepted is True:
            return self.proposal_traj
        elif self.accepted is False:
            return self.initial_traj
        else:
            return None