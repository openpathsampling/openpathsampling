'''
Created on 19.07.2014

@author: Jan-Hendrik Prinz
'''

import numpy as np

import opentis as ops
from opentis.todict import dictable
from opentis.storage.decorators import storable


@dictable
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
        for key, value in kwargs.iteritems():
            setattr(self, key, value)

    def __str__(self):
        # primarily for debugging/interactive use
        mystr = ""
        for key in self.__dict__.keys():
            mystr += str(key) + " = " + str(self.__dict__[key]) + '\n'
        return mystr


@dictable
@storable
class PathMover(object):
    """
    A PathMover is the description of how to generate a new path from an old one.
    
    Notes
    -----
    
    Basically this describes the proposal step for a MC in path space.
    
    We might detach this from the acceptance step?!?!?
    This would mean that a PathMover needs only an old trajectory and gives
    a new one.
    
    For example a ForwardShoot then uses a shooting point selector and runs
    a new trajectory and combine them to get a new one.
    
    After the move has been made, we can retrieve information about the
    move, as well as the new trajectory from the PathMover object
    
    Attributes
    ----------
    engine : DynamicsEngine
        the attached engine used to generate new trajectories
    """

    engine = None

    @property
    def identifier(self):
        if hasattr(self, 'json'):
            return self.json
        else:
            return None

    def __init__(self):
        
        # An ensemble that is at the same time triggering the stopping
        # criterion and the final acceptance and the end.  This is because
        # the goal is to sample trajectories in a specific ensemble. So we
        # want to generate and stop as soon as this cannot be fulfilled
        # anymore. Some of the conditions cannot be checked during runtime,
        # so we have to do that at the end to make sure.

        self.ensemble = ops.FullEnsemble()
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
        After this command additional information can be accessed from this
        object
        '''

        details = MoveDetails(inputs = [trajectory], result = trajectory, success=True)
        path = ops.Sample(trajectory=trajectory, mover=self, ensemble=self.ensemble, details=details)

        return path

    def selection_probability_ratio(self, details=None):
        '''
        Return the proposal probability necessary to correct for an asymmetric proposal.
        
        Notes
        -----
        This is effectively the ratio of proposal probabilities for a mover.
        For symmetric proposal this is one. In the case of e.g. Shooters
        this depends on the used ShootingPointSelector and the start and
        final trajectory.
        
        I am not sure if it makes sense that to define it this way, but for
        Shooters this is, what we need for the acceptance step in addition
        to the check if we have a trajectory of the target ensemble.

        What about Minus Move and PathReversalMove?
        '''
        return 1.0

@dictable
class ShootMover(PathMover):
    '''
    A pathmover that implements a general shooting algorithm that generates
    a sample from a specified ensemble
    '''

    def __init__(self, selector, ensemble):
        super(ShootMover, self).__init__()
        self.selector = selector
        self.ensemble = ensemble
        if PathMover.engine is not None:
            # ask the engine. But when loading this might not exist, so until we store a pointer to the used engine
            # we will create our own stopper or (we say Movers are not creatable but will be stored as stubs
            # so called loaded objects that have the same class name, the same dict, but no functions to call!
            # this will be indicated by removing the @creatable
            self._length_stopper = PathMover.engine.max_length_stopper

    def selection_probability_ratio(self, details):
        '''
        Return the proposal probability for Shooting Moves. These are given
        by the ratio of partition functions
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

        path = ops.Sample(trajectory=details.result, ensemble=self.ensemble, details=details)

        return path
    
    
@dictable
class ForwardShootMover(ShootMover):
    '''
    A pathmover that implements the forward shooting algorithm
    '''
    def _generate(self, details):
        shooting_point = details.start_point.index
        print "Shooting forward from frame %d" % shooting_point
        
        # Run until one of the stoppers is triggered
        partial_trajectory = PathMover.engine.generate(
            details.start_point.snapshot.copy(),
            running = [
                ops.ForwardAppendedTrajectoryEnsemble(
                    self.ensemble, 
                    details.start[0:details.start_point.index]
                ).can_append, 
                self._length_stopper.can_append
            ]
        )

        # DEBUG
        #setattr(details, 'repeated_partial', details.start[0:shooting_point])
        #setattr(details, 'new_partial', partial_trajectory)

        details.final = details.start[0:shooting_point] + partial_trajectory
        details.final_point = ops.ShootingPoint(self.selector, details.final,
                                            shooting_point)
        pass
    
@dictable
class BackwardShootMover(ShootMover):
    '''
    A pathmover that implements the backward shooting algorithm
    '''
    def _generate(self, details):
        print "Shooting backward from frame %d" % details.start_point.index

        # Run until one of the stoppers is triggered
        partial_trajectory = PathMover.engine.generate(
            details.start_point.snapshot.reversed_copy(), 
            running = [
                ops.BackwardPrependedTrajectoryEnsemble(
                    self.ensemble, 
                    details.start[details.start_point.index + 1:]
                ).can_prepend, 
                self._length_stopper.can_prepend
            ]
        )

        # DEBUG
        #setattr(details, 'repeated_partial', details.start[details.start_point.index+1:])
        #setattr(details, 'new_partial', partial_trajectory.reversed)

        details.final = partial_trajectory.reversed + details.start[details.start_point.index + 1:]
        details.final_point = ops.ShootingPoint(self.selector, details.final, partial_trajectory.frames - 1)
        
        pass

@dictable
class MixedMover(PathMover):
    '''
    Defines a mover that picks a over from a set of movers with specific weights.
    
    Notes
    -----
    Channel functions from self.mover to self. Does NOT work yet. Think about a good way to implement this...
    '''
    def __init__(self, movers, weights = None, ensemble = None):
        super(MixedMover, self).__init__()

        self.movers = movers
        if ensemble is not None:
            self.ensemble = ensemble

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

        path = ops.Sample(trajectory=sample.trajectory, ensemble=mover.ensemble, details=sample.details)
        return path

#############################################################
# The following moves still need to be implemented. Check what excactly they do
#############################################################

@dictable
class MinusMove(PathMover):
    def move(self, allpaths, state):
        pass

@dictable
class PathReversal(PathMover):
    def move(self, trajectory, ensemble):
        details = MoveDetails()
        reversed_trajectory = trajectory.reversed()
        details.inputs = [trajectory]
        details.mover = self
        details.final = reversed_trajectory
        details.success = True
        details.acceptance = 1.0
        details.result = reversed_trajectory

        sample = ops.Sample(
            trajectory=details.result,
            mover=self,
            ensemble=ensemble,
            details=details
        )


#############################################################
# The following move should be moved to RETIS and just uses moves. It is not a move itself
#############################################################

@dictable
class ReplicaExchange(PathMover):
    # TODO: Might put the target ensembles into the Mover instance, which means we need lots of mover instances for all ensemble switches
    def move(self, trajectory1, trajectory2, ensemble1, ensemble2):
        success = True # Change to actual check for swapping
        details1 = MoveDetails()
        details2 = MoveDetails()
        details1.inputs = [trajectory1, trajectory2]
        details2.inputs = [trajectory1, trajectory2]
        setattr(details1, 'ensembles', [ensemble1, ensemble2])
        setattr(details2, 'ensembles', [ensemble1, ensemble2])
        details1.mover = self
        details2.mover = self
        details2.final = trajectory1
        details1.final = trajectory2
        if success:
            # Swap
            details1.success = True
            details2.success = True
            details1.acceptance = 1.0
            details2.acceptance = 1.0
            details1.result = trajectory2
            details2.result = trajectory1
        else:
            # No swap
            details1.success = False
            details2.success = False
            details1.acceptance = 0.0
            details2.acceptance = 0.0
            details1.result = trajectory1
            details2.result = trajectory2

        sample1 = ops.Sample(
            trajectory=details1.result,
            mover=self,
            ensemble=ensemble1,
            details=details1
        )
        sample2 = ops.Sample(
            trajectory=details2.result,
            mover=self,
            ensemble=ensemble2,
            details=details2
            )
        return [sample1, sample2]


class PathMoverFactory(object):
    @staticmethod
    def OneWayShootingSet(selector_set, interface_set):
        if type(selector_set) is not list:
            selector_set = [selector_set]*len(interface_set)
        mover_set = [
            MixedMover(movers=[ForwardShootMover(sel, iface), 
                               BackwardShootMover(sel, iface) 
                              ],
                       ensemble=iface
                      )
            for (sel, iface) in zip(selector_set, interface_set)
        ]
        return mover_set

    @staticmethod
    def TwoWayShootingSet():
        pass

    @staticmethod
    def NearestNeighborRepExSet():
        pass
