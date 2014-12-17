'''
Created on 19.07.2014

@author: Jan-Hendrik Prinz
'''

import numpy as np
import random

from shooting import ShootingPoint
from ensemble import ForwardAppendedTrajectoryEnsemble, BackwardPrependedTrajectoryEnsemble
from ensemble import FullEnsemble
from sample import Sample
from wrapper import storable

def make_list_of_pairs(l):
    '''
    Converts input from several possible formats into a list of pairs: used
    to clean input for swap-like moves.

    Allowed input formats: 
    * flat list of length 2N
    * list of pairs
    * None (returns None)

    Anything else will lead to a ValueError or AssertionError
    '''
    if l is None:
        return None

    len_l = len(l) # raises TypeError, avoids everything else

    # based on first element, decide whether this should be a list of lists
    # or a flat list
    try:
        len_l0 = len(l[0])
        list_of_lists = True
    except TypeError:
        list_of_lists = False

    if list_of_lists == True:
        for elem in l:
            assert len(elem)==2, "List of lists: inner list length != 2"
        outlist = l
    else:
        assert len(l) % 2 == 0, "Flattened list: length not divisible by 2"
        outlist = [ [a, b] 
                   for (a, b) in zip(l[slice(0,None,2)], l[slice(1,None,2)])
                  ]
    # Note that one thing we don't check is whether the items are of the
    # same type. That might be worth doing someday; for now, we trust that
    # part to work.
    return outlist


@storable
class MoveDetails(object):
    '''Details of the move as applied to a given replica

    Attributes
    ----------
    replica : integer
        replica ID to which this trial move would apply
    inputs : list of Trajectry
        the Samples which were used as inputs to the move
    final : Tractory
        the Trajectory 
    accepted : bool
        whether the attempted move was accepted
    mover : PathMover
        the PathMover which generated this trial

    Specific move types may have add several other attributes for each
    MoveDetails object. For example, shooting moves will also include
    information about the shooting point selection, etc.

    TODO (or at least to put somewhere):
    rejection_reason : String
        explanation of reasons the path was rejected

    RENAME: inputs=>initial
            final=>trial
            success=>accepted
            accepted=>trial_in_ensemble

    TODO:
    Currently inputs/final/accepted are in terms of Trajectory objects. I
    think it makes more sense for them to be Samples.
    '''

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
    This would mean that a PathMover needs only an old trajectory and gives
    a new one.
    
    For example a ForwardShoot then uses a shooting point selector and runs
    a new trajectory and combine them to get a new one.
    
    After the move has been made, we can retrieve information about the
    move, as well as the new trajectory from the PathMover object
    
    Potential future change: `engine` is not needed for all PathMovers
    (replica exchange, ensemble hopping, path reversal, and moves which
    combine these [state swap] have no need for the engine). Maybe that
    should be moved into only the ensembles that need it? ~~~DWHS

    Also, I agree with the separating trial and acceptance. We might choose
    to use a different acceptance criterion than Metropolis. For example,
    the "waste recycling" approach recently re-discovered by Frenkel (see
    also work by Athenes, Jourdain, and old work by Kalos) might be
    interesting. I think the best way to do this is to keep the acceptance
    in the PathMover, but have it be a ~~~DWHS


    Attributes
    ----------
    engine : DynamicsEngine
        the attached engine used to generate new trajectories

    """

    cls = 'pathmover'
    engine = None

    @property
    def identifier(self):
        if hasattr(self, 'json'):
            return self.json
        else:
            return None

    def __init__(self, replicas='all', ensembles=None):
        
        self.name = self.__class__.__name__

        if type(replicas) is int:
            self.replicas = [replicas]
        else:
            self.replicas = replicas

        if ensembles is not None and type(ensembles) is not list:
            ensembles = [ensembles]
        self.ensembles = ensembles

        self.idx = dict()

    def legal_sample_set(self, globalstate, ensembles=None):
        '''
        This returns all the samples from globalstate which are in both
        self.replicas and the parameter ensembles. If ensembles is None, we
        use self.ensembles. If you want all ensembles allowed, pass
        ensembles='all'.
        '''
        if self.replicas == 'all':
            reps = globalstate.replica_list()
        else:
            reps = self.replicas
        rep_samples = []
        for rep in reps:
            rep_samples.extend(globalstate.all_from_replica(rep))

        if ensembles is None:
            if self.ensembles is None:
                ensembles = 'all'
            else:
                ensembles = self.ensembles

        if ensembles == 'all':
            legal_samples = rep_samples
        else:
            ens_samples = []
            if type(ensembles) is not list:
                ensembles = [ensembles]
            for ens in ensembles:
                ens_samples.extend(globalstate.all_from_ensemble(ens))
            legal_samples = list(set(rep_samples) & set(ens_samples))

        return legal_samples

    def select_sample(self, globalstate, ensembles=None):
        '''
        Returns one of the legal samples given self.replica and the ensemble
        set in ensembles.
        '''
        legal = self.legal_sample_set(globalstate, ensembles)
        return random.choice(legal)

    def move(self, globalstate):
        '''
        Run the generation starting with the initial trajectory specified.

        Parameters
        ----------
        globalstate : GlobalState
            the initial global state
        
        Returns
        -------        
        samples : list of Sample()
            the new samples
        
        Notes
        -----
        After this command additional information can be accessed from this
        object (??? can you explain this, JHP?)
        '''

        return [] # pragma: no cover

    def selection_probability_ratio(self, details=None):
        '''
        Return the proposal probability necessary to correct for an
        asymmetric proposal.
        
        Notes
        -----
        This is effectively the ratio of proposal probabilities for a mover.
        For symmetric proposal this is one. In the case of e.g. Shooters
        this depends on the used ShootingPointSelector and the start and
        final trajectory.
        
        I am not sure if it makes sense that to define it this way, but for
        Shooters this is, what we need for the acceptance step in addition
        to the check if we have a trajectory of
        the target ensemble.

        What about Minus Move and PathReversalMove?
        '''
        return 1.0 # pragma: no cover

class ShootMover(PathMover):
    '''
    A pathmover that implements a general shooting algorithm that generates
    a sample from a specified ensemble 
    '''

    def __init__(self, selector, ensembles=None, replicas='all'):
        super(ShootMover, self).__init__(ensembles=ensembles, replicas=replicas)
        self.selector = selector
        self.length_stopper = PathMover.engine.max_length_stopper
        self.extra_details = ['start', 'start_point', 'final',
                              'final_point']

    def selection_probability_ratio(self, details):
        '''
        Return the proposal probability for Shooting Moves. These are given
        by the ratio of partition functions
        '''
        return details.start_point.sum_bias / details.final_point.sum_bias
    
    def _generate(self, ensemble):
        self.final = self.start
    
    def move(self, globalstate):
        # select a legal sample, use it to determine the trajectory and the
        # ensemble needed for the dynamics
        rep_sample = self.select_sample(globalstate, self.ensembles) 
        trajectory = rep_sample.trajectory
        dynamics_ensemble = rep_sample.ensemble
        replica = rep_sample.replica

        details = MoveDetails()
        details.success = False
        details.inputs = [trajectory]
        details.mover = self
        setattr(details, 'start', trajectory)
        setattr(details, 'start_point', self.selector.pick(details.start) )
        setattr(details, 'final', None)
        setattr(details, 'final_point', None)

        self._generate(details, dynamics_ensemble)


        details.accepted = dynamics_ensemble(details.final)
        details.result = details.start

        if details.accepted:
            rand = np.random.random()
            print 'Proposal probability', self.selection_probability_ratio(details), '/ random :', rand
            if (rand < self.selection_probability_ratio(details)):
                details.success = True
                details.result = details.final

        path = Sample(replica=replica,
                      trajectory=details.result, 
                      ensemble=dynamics_ensemble,
                      details=details)

        return path
    
    
class ForwardShootMover(ShootMover):    
    '''
    A pathmover that implements the forward shooting algorithm
    '''
    def _generate(self, details, ensemble):
        shooting_point = details.start_point.index
        print "Shooting forward from frame %d" % shooting_point
        
        # Run until one of the stoppers is triggered
        partial_trajectory = PathMover.engine.generate(
            details.start_point.snapshot,
            running = [ForwardAppendedTrajectoryEnsemble(
                ensemble, 
                details.start[0:details.start_point.index]).can_append, 
                self.length_stopper.can_append]
             )

        details.final = details.start[0:shooting_point] + partial_trajectory
        details.final_point = ShootingPoint(self.selector, details.final,
                                            shooting_point)
    
class BackwardShootMover(ShootMover):    
    '''
    A pathmover that implements the backward shooting algorithm
    '''
    def _generate(self, details, ensemble):
        print "Shooting backward from frame %d" % details.start_point.index

        # Run until one of the stoppers is triggered
        partial_trajectory = PathMover.engine.generate(
            details.start_point.snapshot.reversed_copy(), 
            running = [BackwardPrependedTrajectoryEnsemble( 
                ensemble, 
                details.start[details.start_point.index + 1:]).can_prepend, 
                self.length_stopper.can_prepend]
        )

        details.final = partial_trajectory.reversed + details.start[details.start_point.index + 1:]
        details.final_point = ShootingPoint(self.selector, details.final, partial_trajectory.frames - 1)
        
        pass

class MixedMover(PathMover):
    '''
    Defines a mover that picks a over from a set of movers with specific
    weights.
    
    Notes
    -----
    Channel functions from self.mover to self. Does NOT work yet. Think
    about a good way to implement this...
    '''
    def __init__(self, movers, ensembles=None, replicas='all', weights = None):
        super(MixedMover, self).__init__(ensembles=ensembles, replicas=replicas)

        self.movers = movers

        if weights is None:
            self.weights = [1.0] * len(movers)
        else:
            self.weights = weights
    
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

        # why do we make a new sample here?
        path = Sample(trajectory=sample.trajectory,
                      ensemble=sample.ensemble, 
                      details=sample.details,
                     replica=sample.replica)
        return path

class ReplicaIDChange(PathMover):
    def __init__(self, new_replicas=None, old_samples=None, 
                 ensembles=None, replicas='all'):
        super(ReplicaIDChange, self).__init__(ensembles, replicas)
        self.new_replicas = new_replicas
        self.old_samples = old_samples

    def move(self, globalstate):
        rep_sample = self.select_sample(globalstate, self.ensembles)
        new_rep = self.new_replicas[rep_sample.replica]
        old_sample = self.old_samples[rep_sample.replica]
        trajectory = rep_sample.trajectory
        
        details = MoveDetails()
        details.inputs = rep_sample.trajectory
        # TODO: details
        dead_sample = Sample(replica=rep_sample.replica,
                             ensemble=old_sample.ensemble,
                             trajectory=old_sample.trajectory
                            )
        new_sample = Sample(replica=new_rep,
                            ensemble=rep_sample.ensemble,
                            trajectory=rep_sample.trajectory
                           )
        return [dead_sample, new_sample]


class EnsembleHopMover(PathMover):
    def __init__(self, bias=None, ensembles=None, replicas='all'):
        ensembles = make_list_of_pairs(ensembles)
        super(EnsembleHopMover, self).__init__(ensembles=ensembles, 
                                               replicas=replicas)
        # TODO: add support for bias: should be a list, one per pair of
        # ensembles -- another version might take a value for each ensemble,
        # and use the ratio; this latter is better for CITIS
        self.bias = bias

    def move(self, globalstate):
        # ensemble hops are in the order [from, to]
        ens_pair = random.choice(self.ensembles)
        ens_from = ens_pair[0]
        ens_to = ens_pair[1]

        rep_sample = self.select_sample(globalstate, ens_from)
        replica = rep_sample.replica
        trajectory = rep_sample.trajectory

        details = MoveDetails()
        details.success = False
        details.inputs = [trajectory]
        details.mover = self
        details.result = trajectory
        setattr(details, 'initial_ensemble', ens_from)
        setattr(details, 'trial_ensemble', ens_to)
        details.success = ens_to(trajectory)
        if details.success == True:
            setattr(details, 'result_ensemble', ens_to)
        else: 
            setattr(details, 'result_ensemble', ens_from)

        path = Sample(trajectory=trajectory,
                      ensemble=details.result_ensemble, 
                      details=details,
                      replica=replica
                     )
        return path



#############################################################
# The following moves still need to be implemented. Check what excactly they do
#############################################################

class MinusMove(PathMover):
    def move(self, allpaths, state):
        pass

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

        sample = Sample(
            trajectory=details.result,
            mover=self,
            ensemble=ensemble,
            details=details
        )


#############################################################
# The following move should be moved to RETIS and just uses moves. It is not a move itself
#############################################################

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

        sample1 = Sample(
            trajectory=details1.result,
            mover=self,
            ensemble=ensemble1,
            details=details1
        )
        sample2 = Sample(
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
            MixedMover([
                ForwardShootMover(sel, ensembles=[iface]), 
                BackwardShootMover(sel, ensembles=[iface])
            ],
            ensembles=[iface])
            for (sel, iface) in zip(selector_set, interface_set)
        ]
        return mover_set

    @staticmethod
    def TwoWayShootingSet():
        pass

    @staticmethod
    def NearestNeighborRepExSet():
        pass
