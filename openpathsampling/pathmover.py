'''
Created on 19.07.2014

@author: Jan-Hendrik Prinz
'''

import numpy as np
import random

import openpathsampling as paths
from openpathsampling.todict import restores_as_stub_object
from sample import Sample, SampleSet
# TODO: switch usage of the next paths.*
from ensemble import Ensemble
from trajectory import Trajectory


import logging
from ops_logging import initialization_logging

logger = logging.getLogger(__name__)
init_log = logging.getLogger('openpathsampling.initialization')

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

@restores_as_stub_object
class MoveDetails(object):
    '''Details of the move as applied to a given replica

    Attributes
    ----------
    replica : integer
        replica ID to which this trial move would apply
    inputs : list of Trajectry
        the Samples which were used as inputs to the move
    trial : Trajectory
        the Trajectory 
    trial_is_in_ensemble : bool
        whether the attempted move created a trajectory in the right
        ensemble
    mover_path : list of PathMover
        the sequence of calls to the PathMover which generated this trial

    Specific move types may have add several other attributes for each
    MoveDetails object. For example, shooting moves will also include
    information about the shooting point selection, etc.

    TODO (or at least to put somewhere):
    rejection_reason : String
        explanation of reasons the path was rejected

    RENAME: inputs=>initial
            accepted=>trial_in_ensemble (probably only in shooting)

    TODO:
    Currently inputs/trial/accepted are in terms of Trajectory objects. I
    think it makes more sense for them to be Samples.
    '''

    def __init__(self, **kwargs):
        self.inputs=None
        self.trial=None
        self.result=None
        self.acceptance_probability=None
        self.accepted=None
        self.mover_path=[]
        for key, value in kwargs:
            setattr(self, key, value)

    def __str__(self):
        # primarily for debugging/interactive use
        mystr = ""
        for key in self.__dict__.keys():
            if not isinstance(self.__dict__[key], paths.Ensemble):
                mystr += str(key) + " = " + str(self.__dict__[key]) + '\n'
        return mystr

    @staticmethod
    def initialization(sample):
        details = MoveDetails()
        details.accepted = True
        details.acceptance_probability = 1.0
        details.mover_path = []
        #details.mover = PathMover()
        #details.mover.name = "Initialization (trajectory)"
        details.inputs = []
        details.trial = sample.trajectory
        details.ensemble = sample.ensemble
        details.result = sample.trajectory
        return details



@restores_as_stub_object
class PathMover(object):
    """
    A PathMover is the description of how to generate a new path from an old
    one.
    
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
    in the PathMover, but have it be a separate class ~~~DWHS


    Attributes
    ----------
    engine : DynamicsEngine
        the attached engine used to generate new trajectories

    """
    engine = None

    def __init__(self, replicas='all', ensembles=None):
        self.name = self.__class__.__name__

        if type(replicas) is int:
            self.replicas = [replicas]
        else:
            self.replicas = replicas

        if ensembles is not None and type(ensembles) is not list:
            ensembles = [ensembles]
        self.ensembles = ensembles

        initialization_logging(logger=init_log, obj=self,
                               entries=['replicas', 'ensembles'])

    def __call__(self, sample_set):
        return sample_set

    def legal_sample_set(self, globalstate, ensembles=None, replicas='all'):
        '''
        This returns all the samples from globalstate which are in both
        self.replicas and the parameter ensembles. If ensembles is None, we
        use self.ensembles. If you want all ensembles allowed, pass
        ensembles='all'.

        TODO: Turn into a filter decorator or
        '''
        if self.replicas == 'all':
            mover_replicas = globalstate.replica_list()
        else:
            mover_replicas = self.replicas
        if replicas == 'all':
            selected_replicas = globalstate.replica_list()
        else:
            selected_replicas = replicas

        reps = list(set(mover_replicas) & set(selected_replicas))
        rep_samples = []
        for rep in reps:
            rep_samples.extend(globalstate.all_from_replica(rep))

        #logger.debug("ensembles = " + str([ensembles]))
        #logger.debug("self.ensembles = " + str(self.ensembles))
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
                #try:
                    #ens_samples.extend(globalstate.all_from_ensemble(ens[0]))
                #except TypeError:
                ens_samples.extend(globalstate.all_from_ensemble(ens))
            legal_samples = list(set(rep_samples) & set(ens_samples))

        return legal_samples

    def select_sample(self, globalstate, ensembles=None, replicas=None):
        '''
        Returns one of the legal samples given self.replica and the ensemble
        set in ensembles.

        TODO: This must be saved somehow (it is actually I think), otherwise
        Samples are not reproducible when applied to a SampleSet!
        '''
        if replicas is None:
            replicas=self.replicas
        logger.debug("replicas: "+str(replicas)+" ensembles: "+repr(ensembles))
        legal = self.legal_sample_set(globalstate, ensembles, replicas)
        for sample in legal:
            logger.debug("legal: (" + str(sample.replica)
                         + "," + str(sample.trajectory)
                         + "," + repr(sample.ensemble)
                         + ")")
        selected = random.choice(legal)
        logger.debug("selected sample: (" + str(selected.replica)
                     + "," + str(selected.trajectory)
                     + "," + repr(selected.ensemble)
                     + ")")
        return selected

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

        return paths.EmptyMovePath() # pragma: no cover

    def selection_probability_ratio(self, details=None):
        '''
        Return the proposal probability necessary to correct for an
        asymmetric proposal.
        
        Notes
        -----
        This is effectively the ratio of proposal probabilities for a mover.
        For symmetric proposal this is one. In the case of e.g. Shooters
        this depends on the used ShootingPointSelector and the start and
        trial trajectory.
        
        I am not sure if it makes sense that to define it this way, but for
        Shooters this is, what we need for the acceptance step in addition
        to the check if we have a trajectory of
        the target ensemble.

        What about Minus Move and PathReversalMove?
        '''
        return 1.0 # pragma: no cover

@restores_as_stub_object
class CollapseMove(PathMover):
    def __init__(self, inner_mover):
        self.inner_mover = inner_mover

    def move(self, globalstate):
        return self.inner_mover.move(globalstate).closed

@restores_as_stub_object
class ShootMover(PathMover):
    '''
    A pathmover that implements a general shooting algorithm that generates
    a sample from a specified ensemble 
    '''

    def __init__(self, selector, ensembles=None, replicas='all'):
        super(ShootMover, self).__init__(ensembles=ensembles,
                                         replicas=replicas
                                        )
        self.selector = selector
        self._length_stopper = PathMover.engine.max_length_stopper
        self._extra_details = ['start', 'start_point', 'trial',
                              'final_point']
        initialization_logging(logger=init_log, obj=self,
                               entries=['selector'])

    def selection_probability_ratio(self, details):
        '''
        Return the proposal probability for Shooting Moves. These are given
        by the ratio of partition functions
        '''
        return details.start_point.sum_bias / details.final_point.sum_bias
    
    def _generate(self, ensemble):
        self.trial = self.start
    
    def move(self, globalstate):
        # select a legal sample, use it to determine the trajectory and the
        # ensemble needed for the dynamics
        rep_sample = self.select_sample(globalstate, self.ensembles)
        trajectory = rep_sample.trajectory
        dynamics_ensemble = rep_sample.ensemble
        replica = rep_sample.replica

        details = MoveDetails()
        details.accepted = False
        details.inputs = [rep_sample]
        details.mover_path.append(self)
        details.mover = self
        setattr(details, 'start', trajectory)
        setattr(details, 'start_point', self.selector.pick(details.start) )
        setattr(details, 'final_point', None)

        self._generate(details, dynamics_ensemble)

        setattr(details, 'trial_is_in_ensemble',
                dynamics_ensemble(details.trial))

        details.result = details.start

        if details.trial_is_in_ensemble:
            rand = np.random.random()
            sel_prob = self.selection_probability_ratio(details)
            logger.info('Proposal probability ' + str(sel_prob)
                        + ' / random : ' + str(rand)
                       )
            if (rand < self.selection_probability_ratio(details)):
                logger.info("Shooting move accepted!")
                details.accepted = True
                details.result = details.trial

        sample = paths.Sample(
            replica=replica, 
            trajectory=details.result, 
            ensemble=dynamics_ensemble,
            details=details
        )

#        new_set = SampleSet(samples=[sample], predecessor=globalstate, accepted=True)
#        new_set = globalstate.apply([sample], accepted = details.accepted, move=self)

        path = paths.SampleMovePath(
                samples=[sample],
                accepted=details.accepted,
                mover=self
        )

        return path
    
    
@restores_as_stub_object
class ForwardShootMover(ShootMover):
    '''
    A pathmover that implements the forward shooting algorithm
    '''
    def _generate(self, details, ensemble):
        shooting_point = details.start_point.index
        shoot_str = "Shooting {sh_dir} from frame {fnum} in [0:{maxt}]"
        logger.info(shoot_str.format(fnum=details.start_point.index,
                                     maxt=len(details.start)-1,
                                     sh_dir="forward"
                                    ))
        
        # Run until one of the stoppers is triggered
        partial_trajectory = PathMover.engine.generate(
            details.start_point.snapshot.copy(),
            running = [
                paths.ForwardAppendedTrajectoryEnsemble(
                    ensemble, 
                    details.start[0:details.start_point.index]
                ).can_append, 
                self._length_stopper.can_append
            ]
        )

        # DEBUG
        #setattr(details, 'repeated_partial', details.start[0:shooting_point])
        #setattr(details, 'new_partial', partial_trajectory)

        details.trial = details.start[0:shooting_point] + partial_trajectory
        details.final_point = paths.ShootingPoint(self.selector, details.trial,
                                            shooting_point)
    
@restores_as_stub_object
class BackwardShootMover(ShootMover):
    '''
    A pathmover that implements the backward shooting algorithm
    '''
    def _generate(self, details, ensemble):
        shoot_str = "Shooting {sh_dir} from frame {fnum} in [0:{maxt}]"
        logger.info(shoot_str.format(fnum=details.start_point.index,
                                     maxt=len(details.start)-1,
                                     sh_dir="backward"
                                    ))

        # Run until one of the stoppers is triggered
        partial_trajectory = PathMover.engine.generate(
            details.start_point.snapshot.reversed_copy(), 
            running = [
                paths.BackwardPrependedTrajectoryEnsemble(
                    ensemble, 
                    details.start[details.start_point.index + 1:]
                ).can_prepend, 
                self._length_stopper.can_prepend
            ]
        )

        # DEBUG
        #setattr(details, 'repeated_partial', details.start[details.start_point.index+1:])
        #setattr(details, 'new_partial', partial_trajectory.reversed)

        details.trial = partial_trajectory.reversed + details.start[details.start_point.index + 1:]
        details.final_point = paths.ShootingPoint(self.selector, details.trial, partial_trajectory.frames - 1)

        pass

@restores_as_stub_object
class RandomChoiceMover(PathMover):
    '''
    Chooses a random mover from its movers list, and runs that move. Returns
    the number of samples the submove return.

    For example, this would be used to select a specific replica exchange
    such that each replica exchange is its own move, and which swap is
    selected at random.

    Attributes
    ----------
    movers : list of PathMover
        the PathMovers to choose from
    weights : list of floats
        the relative weight of each PathMover (does not need to be normalized)
    '''
    def __init__(self, movers, ensembles=None, replicas='all', weights=None, name=None):
        super(RandomChoiceMover, self).__init__(ensembles=ensembles, replicas=replicas)

        if name is not None:
            self.name = name

        self.movers = movers

        if weights is None:
            self.weights = [1.0] * len(movers)
        else:
            self.weights = weights

        initialization_logging(init_log, self,
                               entries=['movers', 'weights'])
    
    def move(self, globalstate):
        rand = np.random.random() * sum(self.weights)
        idx = 0
        prob = self.weights[0]
        while prob <= rand and idx < len(self.weights):
            idx += 1
            prob += self.weights[idx]

        logger_str = "RandomChoiceMover ({name}) selecting mover index {idx} ({mtype})"
        logger.info(logger_str.format(name=self.name, idx=idx, mtype=self.movers[idx].name))

        mover = self.movers[idx]

        path = paths.RandomChoiceMovePath(mover.move(globalstate), mover=self)

        return path

class ConditionalMover(PathMover):
    '''
    An if-then-else structure for PathMovers.

    Returns a SequentialMovePath of the if_move movepath and the then_move
    movepath (if if_move is accepted) or the else_move movepath (if if_move
    is rejected).
    '''
    def __init__(self, if_mover, then_mover, else_mover, ensembles=None,
                 replicas='all'):
        super(ConditionalMover, self).__init__(ensembles=ensembles,
                                               replicas=replicas
                                              )
        self.if_mover = if_mover
        self.then_mover = then_mover
        self.else_mover = else_mover
        initialization_logging(init_log, self,
                               ['if_mover', 'then_mover', 'else_mover'])

    def move(self, globalstate):
        subglobal = globalstate

        ifclause = self.if_mover.move(subglobal)
        samples = ifclause.samples
        subglobal = subglobal.apply_samples(samples)
        
        if ifclause.accepted:
            if self.then_mover is not None:
                resultclause = self.then_mover.move(subglobal)
            else:
                resultclause = paths.EmptyMovePath()
        else:
            if self.else_mover is not None:
                resultclause = self.else_mover.move(subglobal)
            else:
                resultclause = paths.EmptyMovePath()

        return paths.SequentialMovePath([ifclause, resultclause], mover=self)


@restores_as_stub_object
class SequentialMover(PathMover):
    '''
    Performs each of the moves in its movers list. Returns all samples
    generated, in the order of the mover list.

    For example, this would be used to create a move that does a sequence of
    replica exchanges in a given order, regardless of whether the moves
    succeed or fail.
    '''
    def __init__(self, movers, ensembles=None, replicas='all'):
        super(SequentialMover, self).__init__(ensembles=ensembles,
                                              replicas=replicas,
                                             )
        self.movers = movers
        initialization_logging(init_log, self, ['movers'])

    def move(self, globalstate):
        logger.debug("Starting sequential move")

#        subglobal = SampleSet(self.legal_sample_set(globalstate))

        subglobal = globalstate
        movepaths = []

        for mover in self.movers:
            logger.debug("Starting sequential move step "+str(mover))

            # Run the sub mover
            movepath = mover.move(subglobal)
            samples = movepath.samples
            subglobal = subglobal.apply_samples(samples)
            movepaths.append(movepath)

        return paths.SequentialMovePath(movepaths, mover=self)

@restores_as_stub_object
class PartialAcceptanceSequentialMover(SequentialMover):
    '''
    Performs each move in its movers list until complete or until one is not
    accepted. If any move is not accepted, further moves are not attempted,
    but the previous accepted samples remain accepted.

    For example, this would be used to create a bootstrap promotion move,
    which starts with a shooting move, followed by an EnsembleHop/Replica
    promotion ConditionalSequentialMover. Even if the EnsembleHop fails, the
    accepted shooting move should be accepted.
    '''
    def move(self, globalstate):
        logger.debug("==== BEGINNING " + self.name + " ====")
        subglobal = SampleSet(self.legal_sample_set(globalstate))
        movepaths = []
        for mover in self.movers:
            # NOTE: right now, this doesn't quite work correctly if the
            # submovers are also multimovers (e.g., SequentialMovers). We
            # need a way to see whether the move below considered itself
            # accepted; that could mean that submoves of the submove were
            # rejected but the whole submove was accepted, as with
            # SequentialMovers
            # @JHP: does movepath solve the above comment?
            logger.info(str(self.name)
                        + " starting mover index " + str(self.movers.index(mover) )
                        + " (" + mover.name + ")"
                       )
            # Run the sub mover
            movepath = mover.move(subglobal)
            samples = movepath.samples
            subglobal = subglobal.apply_samples(samples)
            movepaths.append(movepath)
            if not movepath.accepted:
                break

        logger.debug("==== FINISHING " + self.name + " ====")
        return paths.PartialAcceptanceSequentialMovePath(movepaths, mover=self)


@restores_as_stub_object
class ConditionalSequentialMover(SequentialMover):
    '''
    Performs each move in its movers list until complete or until one is not
    accepted. If any move in not accepted, all previous samples are updated
    to have set their acceptance to False.

    For example, this would be used to create a minus move, which consists
    of first a replica exchange and then a shooting (extension) move. If the
    replica exchange fails, the move is aborted before doing the dynamics.

    ConditionalSequentialMover only works if there is a *single* active
    sample per replica.
    '''
    def move(self, globalstate):
        logger.debug("Starting conditional sequential move")

#        subglobal = SampleSet(self.legal_sample_set(globalstate))

        subglobal = globalstate
        movepaths = []

        for mover in self.movers:
            logger.debug("Starting sequential move step "+str(mover))

            # Run the sub mover
            movepath = mover.move(subglobal)
            samples = movepath.samples
            subglobal = subglobal.apply_samples(samples)
            movepaths.append(movepath)

            if not movepath.accepted:
                break

        return paths.ConditionalSequentialMovePath(movepaths, mover=self)


@restores_as_stub_object
class RestrictToLastSampleMover(PathMover):
    def __init__(self, inner_mover):
        self.inner_mover = inner_mover

    def move(self, globalstate):
        movepath = self.inner_mover.move(globalstate)
        return paths.KeepLastSampleMovePath(movepath, mover=self)

@restores_as_stub_object
class ReplicaIDChangeMover(PathMover): 
    """
    Changes the replica ID for a path.
    """
    def __init__(self, replica_pairs, ensembles=None, replicas='all'):
        self.replica_pairs = make_list_of_pairs(replica_pairs)
        super(ReplicaIDChangeMover, self).__init__(ensembles=ensembles,
                                                   replicas=replicas
                                                  )
        self._extra_details = ['rep_from', 'rep_to']
        initialization_logging(logger=init_log, obj=self, 
                               entries=['replica_pairs'])

    def move(self, globalstate):
        legal_from_rep = [rep[0] for rep in self.replica_pairs]
        rep_sample = self.select_sample(globalstate,
                                        ensembles=self.ensembles, 
                                        replicas=legal_from_rep)

        legal_pairs = [pair for pair in self.replica_pairs
                       if pair[0]==rep_sample.replica]
        mypair = random.choice(legal_pairs)

        details = MoveDetails()
        details.inputs = [rep_sample]
        details.trial = rep_sample.trajectory
        details.result = rep_sample.trajectory
        details.accepted = True
        details.acceptance_probability = 1.0
        setattr(details, 'rep_from', mypair[0])
        setattr(details, 'rep_to', mypair[1])

        logger.info("Creating new sample from replica ID " + str(mypair[0])
                    + " and putting it in replica ID " + str(mypair[1]))

        # note: currently this clones into a new replica ID. We might later
        # want to kill the old replica ID (and possibly rename this mover).
        new_sample = paths.Sample(
            replica=details.rep_to, 
            ensemble=rep_sample.ensemble,
            trajectory=rep_sample.trajectory,
            details=details
        )

        return paths.SampleMovePath(
            [new_sample],
            mover=self,
            accepted=details.accepted
        )

@restores_as_stub_object
class EnsembleHopMover(PathMover):
    def __init__(self, bias=None, ensembles=None, replicas='all'):
        # TODO: maybe allow a version of this with a single ensemble and ANY
        # ensemble can hop to that? messy to code; maybe same idea under
        # another name
        ensembles = make_list_of_pairs(ensembles)
        super(EnsembleHopMover, self).__init__(ensembles=ensembles, 
                                               replicas=replicas
                                              )
        # TODO: add support for bias: should be a list, one per pair of
        # ensembles -- another version might take a value for each ensemble,
        # and use the ratio; this latter is better for CITIS
        self.bias = bias
        initialization_logging(logger=init_log, obj=self,
                               entries=['bias'])

    def select_ensemble_pair(self, globalstate):
        # ensemble hops are in the order [from, to]
        initial_ensembles = [pair[0] for pair in self.ensembles]
        logger.debug("initial_ensembles: " + str(initial_ensembles))
        legal_ensembles = [
            s.ensemble
            for s in self.legal_sample_set(globalstate, initial_ensembles)
        ]
        logger.debug("globalstate ensembles" + 
                     str([s.ensemble for s in globalstate]))
        logger.debug("self.ensembles: " + str(self.ensembles))
        logger.debug("Legal Ensembles: " + str(legal_ensembles))
        legal_pairs = [pair for pair in self.ensembles 
                       if pair[0] in legal_ensembles]
        logger.debug("Legal pairs: " + str(legal_pairs))
        ens_pair = random.choice(legal_pairs)
        return ens_pair

    def move(self, globalstate):
        ens_pair = self.select_ensemble_pair(globalstate)
        ens_from = ens_pair[0]
        ens_to = ens_pair[1]

        rep_sample = self.select_sample(globalstate, ens_from)
        logger.debug("Selected sample: " + repr(rep_sample))
        replica = rep_sample.replica

        logger.info("Attempting ensemble hop from {e1} to {e2} replica ID {rid}".format(
            e1=repr(ens_from), e2=repr(ens_to), rid=repr(replica)))

        trajectory = rep_sample.trajectory
        logger.debug("  selected replica: " + str(replica))
        logger.debug("  initial ensemble: " + repr(rep_sample.ensemble))

        details = MoveDetails()
        details.inputs = [rep_sample]
        details.result = trajectory
        setattr(details, 'initial_ensemble', ens_from)
        setattr(details, 'trial_ensemble', ens_to)
        details.accepted = ens_to(trajectory)
        logger.info("Hop starts from legal ensemble: "+str(ens_from(trajectory)))
        logger.info("Hop ends in legal ensemble: "+str(ens_to(trajectory)))

        if details.accepted == True:
            setattr(details, 'result_ensemble', ens_to)
        else:
            setattr(details, 'result_ensemble', ens_from)

        sample = paths.Sample(
            replica=replica,
            trajectory=trajectory,
            ensemble=details.result_ensemble,
            details=details
        )

        path = paths.SampleMovePath(
            [sample],
            mover=self,
            accepted=details.accepted
        )

        return path


@restores_as_stub_object
class ForceEnsembleChangeMover(EnsembleHopMover):
    '''
    Force an ensemble change in the sample.

    This should only be used as part of other moves, since this can create
    samples which are not valid.
    '''
    def __init__(self, ensembles=None, replicas='all'):
        # no bias allowed
        super(ForceEnsembleChangeMover, self).__init__(ensembles=ensembles,
                                                       replicas=replicas)

    def move(self, globalstate):
        ens_pair = self.select_ensemble_pair(globalstate)
        ens_from = ens_pair[0]
        ens_to = ens_pair[1]
        rep_sample = self.select_sample(globalstate, ens_from)
        logger.debug("Selected sample: " + repr(rep_sample))

        replica = rep_sample.replica
        trajectory = rep_sample.trajectory

        details = MoveDetails()
        details.accepted = True
        details.inputs = [rep_sample]
        details.mover_path.append(self)
        details.result = trajectory
        setattr(details, 'initial_ensemble', ens_from)
        setattr(details, 'trial_ensemble', ens_to)
        setattr(details, 'result_ensemble', ens_to)

        sample = paths.Sample(
            trajectory=trajectory,
            ensemble=details.result_ensemble,
            details=details,
            replica=replica
        )

        path = paths.SampleMovePath( [sample], mover=self, accepted=details.accepted)

        return path


@restores_as_stub_object
class RandomSubtrajectorySelectMover(PathMover):
    '''
    Samples a random subtrajectory satisfying the given subensemble.

    If there are no subtrajectories which satisfy the subensemble, this
    returns the zero-length trajectory.
    '''
    def __init__(self, subensemble, n_l=None, ensembles=None, replicas='all'):
        super(RandomSubtrajectorySelectMover, self).__init__(
            ensembles=ensembles, replicas=replicas
        )
        self.n_l=n_l
        self._subensemble = subensemble

    def _choose(self, trajectory_list):
        return random.choice(trajectory_list)

    def move(self, globalstate):
        rep_sample = self.select_sample(globalstate)
        trajectory = rep_sample.trajectory
        replica = rep_sample.replica
        logger.debug("Working with replica " + str(replica) + " (" + str(trajectory) + ")")

        details = MoveDetails()
        details.inputs = [rep_sample]
        details.mover_path.append(self)

        subtrajs = self._subensemble.split(trajectory)
        logger.debug("Found "+str(len(subtrajs))+" subtrajectories.")

#        if self.n_l is None:
#            length_req = lambda x: x > 0
#        else:
#            length_req = lambda x: x==self.n_l

#        if length_req(len(subtrajs)):

        if (self.n_l is None and len(subtrajs) > 0) or \
            (self.n_l is not None and len(subtrajs) == self.n_l):
            subtraj = self._choose(subtrajs)
        else:
            # return zero-length trajectory otherwise
            subtraj = Trajectory([])

        details.trial = subtraj
        details.accepted = True
        details.result = subtraj
        details.acceptance_probability = 1.0
        
        sample = Sample(
            replica=replica,
            trajectory=details.result,
            ensemble=self._subensemble,
            details=details
        )
        path = paths.SampleMovePath( [sample], mover=self, accepted=details.accepted)

        return path

@restores_as_stub_object
class FirstSubtrajectorySelectMover(RandomSubtrajectorySelectMover):
    '''
    Samples the first subtrajectory satifying the given subensemble.

    If there are no subtrajectories which satisfy the ensemble, this returns
    the zero-length trajectory.
    '''
    def _choose(self, trajectory_list):
        return trajectory_list[0]

@restores_as_stub_object
class FinalSubtrajectorySelectMover(RandomSubtrajectorySelectMover):
    '''
    Samples the final subtrajectory satifying the given subensemble.

    If there are no subtrajectories which satisfy the ensemble, this returns
    the zero-length trajectory.
    '''
    def _choose(self, trajectory_list):
        return trajectory_list[-1]

@restores_as_stub_object
class PathReversalMover(PathMover):
    def move(self, globalstate):
        rep_sample = self.select_sample(globalstate, self.ensembles)
        trajectory = rep_sample.trajectory
        ensemble = rep_sample.ensemble
        replica = rep_sample.replica

        details = MoveDetails()
        details.inputs = [rep_sample]
        details.mover_path.append(self)

        reversed_trajectory = trajectory.reversed
        details.trial = reversed_trajectory

        details.accepted = ensemble(reversed_trajectory)
        logger.info("PathReversal move accepted: "+str(details.accepted))
        if details.accepted == True:
            details.acceptance_probability = 1.0
            details.result = reversed_trajectory
        else:
            details.acceptance_probability = 0.0
            details.result = trajectory

        sample = Sample(
            replica=replica,
            trajectory=details.result,
            ensemble=ensemble,
            details=details
        )

        path = paths.SampleMovePath( [sample], mover=self, accepted=details.accepted)

        return path

@restores_as_stub_object
class ReplicaExchangeMover(PathMover):
    def __init__(self, bias=None, ensembles=None, replicas='all'):
        if replicas=='all' and ensembles is None:
            # replicas MUST be a non-empty list of pairs
            raise ValueError("Specify replicas for ReplicaExchangeMover")
        if replicas != 'all':
            replicas = make_list_of_pairs(replicas)
        ensembles = make_list_of_pairs(ensembles)
        # either replicas or ensembles must be a list of pairs; more
        # complicated filtering can be done with a wrapper class
        super(ReplicaExchangeMover, self).__init__(ensembles=ensembles, 
                                                   replicas=replicas)
        # TODO: add support for bias; cf EnsembleHopMover
        self.bias = bias
        initialization_logging(logger=init_log, obj=self,
                               entries=['bias'])


    def move(self, globalstate):
        if self.ensembles is not None:
            [ens1, ens2] = random.choice(self.ensembles)
            s1 = self.select_sample(globalstate, ens1)
            s2 = self.select_sample(globalstate, ens2)
        else:
            [rep1, rep2] = random.choice(self.replicas)
            s1 = globalstate[rep1]
            s2 = globalstate[rep2]
        
        # convert sample to the language used here before
        trajectory1 = s1.trajectory
        trajectory2 = s2.trajectory
        ensemble1 = s1.ensemble
        ensemble2 = s2.ensemble
        replica1 = s1.replica
        replica2 = s2.replica

        from1to2 = ensemble2(trajectory1)
        logger.debug("trajectory " + repr(trajectory1) +
                     " into ensemble " + repr(ensemble2) +
                     " : " + str(from1to2))
        from2to1 = ensemble1(trajectory2)
        logger.debug("trajectory " + repr(trajectory2) +
                     " into ensemble " + repr(ensemble1) +
                     " : " + str(from2to1))
        allowed = from1to2 and from2to1
        details1 = MoveDetails()
        details2 = MoveDetails()
        details1.inputs = [s1, s2]
        details2.inputs = [s2, s1]
        setattr(details1, 'ensembles', [ensemble1, ensemble2])
        setattr(details2, 'ensembles', [ensemble1, ensemble2])
        details1.mover_path.append(self)
        details2.mover_path.append(self)
        details2.trial = trajectory1
        details1.trial = trajectory2
        if allowed:
            # Swap
            details1.accepted = True
            details2.accepted = True
            details1.acceptance_probability = 1.0
            details2.acceptance_probability = 1.0
            details1.result = trajectory2
            details2.result = trajectory1
        else:
            # No swap
            details1.accepted = False
            details2.accepted = False
            details1.acceptance_probability = 0.0
            details2.acceptance_probability = 0.0
            details1.result = trajectory1
            details2.result = trajectory2

        sample1 = paths.Sample(
            replica=replica1,
            trajectory=details1.result,
            ensemble=ensemble1,
            details=details1
        )
        sample2 = paths.Sample(
            replica=replica2,
            trajectory=details2.result,
            ensemble=ensemble2,
            details=details2
            )

        path = paths.SampleMovePath([sample1, sample2], mover=self, accepted=details1.accepted)

        return path


class FilterByReplica(PathMover):
    def __init__(self, mover, replicas):
        if type(replicas) is not list:
            replicas = [replicas]
        self.replicas = replicas
        self.mover = mover
        # TODO: clean this up
        pass

    def move(self, globalstate):
        filtered_gs = SampleSet(
            [s for s in globalstate if s.replica in self.replicas]
        )
        return self.mover.move(filtered_gs)

class OneWayShootingMover(RandomChoiceMover):
    '''
    OneWayShootingMover is a special case of a RandomChoiceMover which
    combines gives a 50/50 chance of selecting either a ForwardShootMover or
    a BackwardShootMover. Both submovers use the same shooting point
    selector, and both apply to the same ensembles and replicas.

    Attributes
    ----------
    sel : ShootingPointSelector
        The shooting point selection scheme
    ensembles : list of Ensemble or None
        valid ensembles; None implies all ensembles are allowed (no
        restriction)
    replicas : list or 'all'
        valid replicas
    '''
    def __init__(self, sel, ensembles=None, replicas='all'):
        movers = [
            ForwardShootMover(sel, ensembles),
            BackwardShootMover(sel, ensembles)
        ]
        super(OneWayShootingMover, self).__init__(
            movers=movers, ensembles=ensembles, replicas=replicas
        )

@restores_as_stub_object
class MinusMover(ConditionalSequentialMover):
    '''
    Instance of a MinusMover.

    The minus move combines a replica exchange with path extension to swap
    paths between the innermost regular TIS interface ensemble and the minus
    interface ensemble. This is particularly useful for improving sampling
    of path space.
    '''
    def __init__(self, minus_ensemble, innermost_ensemble, 
                 ensembles=None, replicas='all'):
        segment = minus_ensemble._segment_ensemble
        subtrajectory_selector = RandomChoiceMover([
            FirstSubtrajectorySelectMover(subensemble=segment,
                                          n_l=minus_ensemble.n_l,
                                          ensembles=[minus_ensemble]
                                         ),
            FinalSubtrajectorySelectMover(subensemble=segment, 
                                          n_l=minus_ensemble.n_l,
                                          ensembles=[minus_ensemble]
                                         ),
        ])
        subtrajectory_selector.name = "MinusSubtrajectoryChooser"

        repex = ReplicaExchangeMover(ensembles=[[segment, innermost_ensemble]])

        force_to_minus = ForceEnsembleChangeMover(
            ensembles=[[segment, minus_ensemble]]
        )

        extension_mover = RandomChoiceMover([
            ForwardShootMover(paths.FinalFrameSelector(), minus_ensemble),
            BackwardShootMover(paths.FirstFrameSelector(), minus_ensemble)
        ])
        extension_mover.name = "MinusExtensionDirectionChooser"

        movers = [
            subtrajectory_selector,
            repex,
            force_to_minus,
            extension_mover
        ]

        self.minus_ensemble = minus_ensemble
        self.innermost_ensemble = innermost_ensemble
        initialization_logging(init_log, self, ['minus_ensemble',
                                                'innermost_ensemble'])
        super(MinusMover, self).__init__(movers=movers,
                                         ensembles=ensembles,
                                         replicas=replicas)

        return

    def move(self, globalstate):
        return super(MinusMover, self).move(globalstate).closed

@restores_as_stub_object
class MultipleSetMinusMover(RandomChoiceMover):
    pass

@restores_as_stub_object
def NeighborEnsembleReplicaExchange(ensemble_list):
    movers = [
        ReplicaExchangeMover(ensembles=[[ensemble_list[i], ensemble_list[i+1]]])
        for i in range(len(ensemble_list)-1)
    ]
    return movers

def PathReversalSet(l):
    if isinstance(l[0], paths.Ensemble):
        return [PathReversalMover(ensembles=[item]) for item in l]
    else:
        return [PathReversalMover(replicas=[item]) for item in l]


class PathMoverFactory(object):
    @staticmethod
    def OneWayShootingSet(selector_set, interface_set):
        if type(selector_set) is not list:
            selector_set = [selector_set]*len(interface_set)
        mover_set = [
            OneWayShootingMover(sel=sel, ensembles=[iface])
            for (sel, iface) in zip(selector_set, interface_set)
        ]
        return mover_set

    @staticmethod
    def TwoWayShootingSet():
        pass

    @staticmethod
    def NearestNeighborRepExSet():
        pass

