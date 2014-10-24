'''
Created on 19.07.2014

@author: Jan-Hendrik Prinz
'''

import numpy as np

from shooting import ShootingPoint
from ensemble import ForwardAppendedTrajectoryEnsemble, BackwardPrependedTrajectoryEnsemble, LengthEnsemble
from storage.wrapper import storable
from ensemble import FullEnsemble
from trajectory import Sample


@storable
class MoveDetails(object):

    def __init__(self, **kwargs):
        self.inputs=None
        self.final=None
        self.result=None
        self.acceptance=None
        self.success=None
        self.accepted=None
        self.mover=None
        for key, value in kwargs:
            setattr(self, key, value)

class PathMover(object):
    """
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
    """

    cls = 'pathmover'
    simulator = None

    @property
    def identifier(self):
        if hasattr(self, 'json'):
            return self.json
        else:
            return None

    def __init__(self):
        
        # An ensemble that is at the same time triggering the stopping criterion and the final acceptance and the end.
        # This is because the goal is to sample trajectories in a specific ensemble. So we want to generate and stop
        # as soon as this cannot be fulfilled anymore. Some of the conditions cannot be checked during runtime, so we
        # have to do that at the end to make sure.

        self.ensemble = FullEnsemble()
        self.name = self.__class__.__name__

        self.idx = dict()

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

        details = MoveDetails(inputs = [trajectory], result = trajectory, success=True)
        path = Sample(trajectory=trajectory, mover=self, ensemble=self.ensemble, details=details)

        return path

    def selection_probability_ratio(self, details=None):
        '''
        Return the proposal probability necessary to correct for an asymmetric proposal.
        
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
        self.length_stopper = PathMover.simulator.max_length_stopper

    def selection_probability_ratio(self, details):
        '''
        Return the proposal probability for Shooting Moves. These are given by the ratio of partition functions
        '''
        return details.start_point.sum_bias / details.final_point.sum_bias
    
    def _generate(self):
        self.final = self.start
    
    def move(self, trajectory):
        details = MoveDetails()
        details.success = False
        details.inputs = [trajectory]
        details.mover = self
        setattr(details, 'start', trajectory)
        setattr(details, 'start_point', self.selector.pick(details.start) )
        setattr(details, 'final', None)
        setattr(details, 'final_point', None)

        self._generate(details)


        details.accepted = self.ensemble(details.final)
        details.result = details.start

        if details.accepted:
            rand = np.random.random()
            print 'Proposal probability', self.selection_probability_ratio(details), '/ random :', rand
            if (rand < self.selection_probability_ratio(details)):
                details.success = True
                details.result = details.final

        path = Sample(trajectory=details.result, mover=self, ensemble=self.ensemble, details=details)

        print path
        return path
    
    
class ForwardShootMover(ShootMover):    
    '''
    A pathmover that implements the forward shooting algorithm
    '''
    def _generate(self, details):
        print "Shooting forward from frame %d" % details.start_point.index
        
        # Run until one of the stoppers is triggered
        partial_trajectory = PathMover.simulator.generate(
                                     details.start_point.snapshot,
                                     running = [ForwardAppendedTrajectoryEnsemble(
                                                      self.ensemble,
                                                      details.start[0:details.start_point.index]
                                                ).forward, self.length_stopper.forward]
                                     )

        details.final = details.start[0:details.start_point.index] + partial_trajectory
        details.final_point = ShootingPoint(self.selector, details.final, details.start_point.index)

        pass
    
class BackwardShootMover(ShootMover):    
    '''
    A pathmover that implements the backward shooting algorithm
    '''
    def _generate(self, details):
        print "Shooting backward from frame %d" % details.start_point.index

        # Run until one of the stoppers is triggered
        partial_trajectory = PathMover.simulator.generate(
                                     details.start_point.snapshot.reversed_copy(),
                                     running = [BackwardPrependedTrajectoryEnsemble(
                                                     self.ensemble,
                                                     details.start[details.start_point.index + 1:]
                                                ).backward, self.length_stopper.backward]
                                     )

        details.final = partial_trajectory.reversed + details.start[details.start_point.index + 1:]
        details.final_point = ShootingPoint(self.selector, details.final, partial_trajectory.frames - 1)
        
        pass

class MixedMover(PathMover):
    '''
    Defines a mover that picks a over from a set of movers with specific weights.
    
    Notes
    -----
    Channel functions from self.mover to self. Does NOT work yet. Think about a good way to implement this...
    '''
    def __init__(self, movers, weights = None):
        super(MixedMover, self).__init__()

        self.movers = movers

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

        mover = self.movers[idx]

        sample = mover.move(trajectory)
        setattr(sample.details, 'mover_idx', idx)

        path = Sample(trajectory=sample.trajectory, mover=self, ensemble=mover.ensemble, details=sample.details)
        return path

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