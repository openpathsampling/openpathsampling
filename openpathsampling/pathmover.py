'''
Created on 19.07.2014

@author: Jan-Hendrik Prinz
'''

import numpy as np
import random

import openpathsampling as paths
from openpathsampling.todict import OPSNamed, OPSObject

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

def keep_selected_samples(func):
    def wrapper(self, *args, **kwargs):
        if 'keep_samples' in kwargs:
            keep_samples = kwargs['keep_samples']
            del kwargs['keep_samples']
            movepath = func(self, *args, **kwargs)
            return paths.FilterSamplesPathMoveChange(movepath, selected_samples=keep_samples)
        else:
            movepath = func(self, *args, **kwargs)
            return movepath

    return wrapper


class PathMover(OPSNamed):
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

    def __init__(self,  ensembles=None):
        super(PathMover, self).__init__()
#        self.name = self.__class__.__name__

        if ensembles is not None and type(ensembles) is not list:
            ensembles = [ensembles]
        self.ensembles = ensembles

        initialization_logging(logger=init_log, obj=self,
                               entries=['ensembles'])

    @property
    def submovers(self):
        return []

    def __iter__(self):
        yield self
        for submove in self.submovers:
            for change in submove:
                yield change

    def __getitem__(self, item):
        if type(item) is int:
            return self.submovers[item]

    def __reversed__(self):
        for submove in self.submovers:
            for change in reversed(submove):
                yield change

        yield self

    def __len__(self):
        if self._len is None:
            self._len = len(list(iter(self)))

        return self._len

    def key(self, change):
        tree = self.keytree()
        return [leave for leave in tree if leave[1] is change ][0][0]

    def _check_head_node(self, items):
        if isinstance(items[0], paths.PathMover):
            # a subtree of pathmovers
            if self.mover is items[0]:
                #print 'found head'
                # found current head node, check, if children match in order
                left = 0
                submovers = [ch.mover for ch in self.submovers]
                subvalues = items[1]
                if type(subvalues) is not list:
                    subvalues = [subvalues]

                for sub in zip(subvalues[0::2], subvalues[1::2]):
                    if left >= len(self.submovers):
                        # no more submoves to match
                        return False
                    if sub is None:
                        # None is a placeholder so move token +1
                        left = left + 1
                    if type(sub) is dict:
                        if sub[0] is None:
                            while left < len(self.submovers):
                                if not [self.submovers[left].mover, [sub[1]]] in self.submovers[left]:
                                    left = left + 1
                                else:
                                    left = left + 1
                                    break

                            if left == len(self.submovers):
                                return False
                        elif sub[0] not in submovers[left:]:
                            #print 'missing sub', sub.keys()[0], 'in', submovers[left:]
                            return False
                        else:
                            idx = submovers.index(sub[0])
                            left = idx + 1
                            if not [sub[0], [sub[1]]] in self.submovers[idx]:
                                #print 'try', {sub.keys()[0] : sub.values()[0]}
                                return False

                    elif isinstance(sub, paths.PathMover):
                        if sub not in submovers[left:]:
                            return False
                        idx = submovers.index(sub)
                        left = idx + 1

                return True

        elif items[0] is None or len(items) == 0:
            # means empty tree and since nothing is in every tree return true
            return True

    def __contains__(self, item):
        """
        Check if a pathmover, pathmovechange or a tree is in self

        A node is either None or a PathMover

        1. Subchanges are given using a dict { parent : child }
        2. Several submoves are given in a list. [child1, child2]
        3. A single submove can be given as a list of length 1 or a single mover.
        4. None is a wildcat and matches everything

        Examples
        --------
        >>> tree1 = {mover1 : mover2}
        >>> tree2 = {mover1 : [mover2, mover3]}
        >>> tree3 = {mover1 : [mover2, {mover4 : [mover5]}] }
        >>> tree4 = {}

        Notes
        -----
        TODO: Add other types of nodes. e.g. explicit PathMoveChange,
        Boolean for .accepted

        Parameters
        ----------
        item : PathMover, PathMoveChange, PathMoveTree

        """
        if isinstance(item, paths.PathMover):
            return item in self.map_post_order(lambda x : x.mover)
        elif type(item) is list:
            if self._check_head_node(item):
                return True

            # Disable checking for submoves for now

            # the head node did not fit so continue trying subnodes
#            for sub in self.submovers:
#                if item in sub:
#                    return True

            return False

        else:
            raise ValueError('Only PathMovers or PathMoveChanges can be tested.')

    def tree(self):
        return {self : [ ch.tree() for ch in self.submovers] }

    def movetree(self):
        return {self.mover : [ ch.movetree() for ch in self.submovers] }

    def keytree(self, movepath=None):

        if movepath is None:
            movepath = [self.mover]

        result = list()
        result.append( ( movepath, self ) )
        mp = []
        for sub in self.submovers:
            subtree = sub.keytree()
            result.extend([ ( movepath + [mp + m[0]], m[1] ) for m in subtree ])
#            print subtree[-1][0]
            mp = mp + subtree[-1][0]


        return result

    def map_tree(self, fnc, **kwargs):
        """
        Apply a function to each node and return the tree

        Parameters
        ----------
        fnc : function(pathmover, args, kwargs)
            the function run at each pathmover node. It is given the node
            and the optional (fixed) parameters
        kwargs : named arguments
            optional arguments added to the function

        Returns
        -------
        tree (fnc(node, **kwargs))
            nested list of the results of the map
        """

        if len(self.submovers) > 1:
            return { fnc(self, **kwargs) : [node.map_tree(fnc, **kwargs) for node in self.submovers]}
        elif len(self.submovers) == 1:
            return { fnc(self, **kwargs) : self.submovers[0].map_tree(fnc, **kwargs)}
        else:
            return fnc(self, **kwargs)

    def map_post_order(self, fnc, **kwargs):
        """
        Traverse the tree of pathmovers in post-order applying a function

        This maps the underlying tree of pathmovers and applies the
        given function at each node returning a list of the results. Post-order
        will result in the order in which samples are generated. That means
        that submoves are called first BEFORE the node itself is evaluated.

        Parameters
        ----------
        fnc : function(pathmover, args, kwargs)
            the function run at each pathmover node. It is given the node
            and the optional (fixed) parameters
        kwargs : named arguments
            optional arguments added to the function

        Returns
        -------
        list (fnc(node, **kwargs))
            flattened list of the results of the map

        Notes
        -----
        This uses the same order as `reversed()`

        See also
        --------
        map_pre_order, map_post_order, level_pre_order, level_post_order
        """
        return [ fnc(node, **kwargs) for node in reversed(self) ]

    def level_post_order(self, fnc, level=0, **kwargs):
        """
        Traverse the tree of pathmovers in post-order applying a function

        This maps the underlying tree of pathmovers and applies the
        given function at each node returning a list of the results. Post-order
        will result in the order in which samples are generated. That means
        that submoves are called first BEFORE the node itself is evaluated.

        Parameters
        ----------
        fnc : function(pathmover, args, kwargs)
            the function run at each pathmover node. It is given the node
            and the optional parameters
        level : int
            the initial level
        kwargs : named arguments
            optional arguments added to the function

        Returns
        -------
        list of tuple(level, func(node, **kwargs))
            flattened list of tuples of results of the map. First part of
            the tuple is the level, second part is the function result.

        See also
        --------
        map_pre_order, map_post_order, level_pre_order, level_post_order
        """

        output = list()
        for mp in self.submovers:
            output.extend(mp.level_post_order(fnc, level + 1, **kwargs))
        output.append((level, fnc(self, **kwargs)))

        return output

    def map_pre_order(self, fnc, **kwargs):
        """
        Traverse the tree of pathmovers in pre-order applying a function

        This maps the underlying tree of pathmovers and applies the
        given function at each node returning a list of the results. Pre-order
        means that submoves are called AFTER the node itself is evaluated.

        Parameters
        ----------
        fnc : function(pathmover, args, kwargs)
            the function run at each pathmover node. It is given the node
            and the optional parameters
        kwargs : named arguments
            optional arguments added to the function

        Returns
        -------
        list (fnc(node, **kwargs))
            flattened list of the results of the map

        Notes
        -----
        This uses the same order as `iter()`

        See also
        --------
        map_pre_order, map_post_order, level_pre_order, level_post_order
        """
        return [ fnc(node, **kwargs) for node in iter(self) ]

    def level_pre_order(self, fnc, level=0, **kwargs):
        """
        Traverse the tree of pathmovers in pre-order applying a function

        This maps the underlying tree of pathmovers and applies the
        given function at each node returning a list of the results. Pre-order
        means that submoves are called AFTER the node itself is evaluated.

        Parameters
        ----------
        fnc : function(pathmover, args, kwargs)
            the function run at each pathmover node. It is given the node
            and the optional parameters
        level : int
            the initial level
        kwargs : named arguments
            optional arguments added to the function

        Returns
        -------
        list of tuple(level, fnc(node, **kwargs))
            flattened list of tuples of results of the map. First part of
            the tuple is the level, second part is the function result.


        See also
        --------
        map_pre_order, map_post_order, level_pre_order, level_post_order
        """

        output = list()
        output.append((level, fnc(self, **kwargs)))

        for mp in self.submovers:
            output.extend(mp.level_pre_order(fnc, level + 1, **kwargs))

        return output

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
        mover_replicas = globalstate.replica_list()

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
            replicas='all'

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

    @keep_selected_samples
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

        return paths.EmptyPathMoveChange() # pragma: no cover

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

    def __str__(self):
        if self.name == self.__class__.__name__:
            return self.__repr__()
        else:
            return self.name


class CollapseMove(PathMover):
    def __init__(self, inner_mover):
        self.inner_mover = inner_mover

    def move(self, globalstate):
        return self.inner_mover.move(globalstate).closed

class ShootMover(PathMover):
    '''
    A pathmover that implements a general shooting algorithm that generates
    a sample from a specified ensemble 
    '''

    def __init__(self, selector, ensembles=None, replicas='all'):
        super(ShootMover, self).__init__(ensembles=ensembles)
        self.selector = selector
        if hasattr(PathMover, 'engine') and hasattr(PathMover.engine, 'max_length_stopper'):
            self._length_stopper = PathMover.engine.max_length_stopper
        else:
            self._length_stopper = paths.FullEnsemble().can_append

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
    
    def _generate(self, details, ensemble):
        self.trial = self.start

    @keep_selected_samples
    def move(self, globalstate):
        # select a legal sample, use it to determine the trajectory and the
        # ensemble needed for the dynamics
        rep_sample = self.select_sample(globalstate, self.ensembles)
        trajectory = rep_sample.trajectory
        dynamics_ensemble = rep_sample.ensemble
        replica = rep_sample.replica

        sample_details = SampleDetails()
        setattr(sample_details, 'start_point', self.selector.pick(trajectory) )
        sample_details.start = trajectory


        self._generate(sample_details, dynamics_ensemble)

        logger.info("Trial trajectory: " +
                    dynamics_ensemble.trajectory_summary_str(sample_details.trial))
        valid = dynamics_ensemble(sample_details.trial)
        accepted = False

        sel_prob = self.selection_probability_ratio(sample_details)
        sample_details.selection_probability = sel_prob

        if valid:
            rand = np.random.random()
            logger.info('Proposal probability ' + str(sel_prob)
                        + ' / random : ' + str(rand)
                       )

            if (rand < sel_prob):
                logger.info("Shooting move accepted!")
                accepted = True

            sample_details.acceptance_probability = sel_prob
        else:
            sample_details.acceptance_probability = 0.0



        trial = paths.Sample(
            replica=replica,
            trajectory=sample_details.trial,
            ensemble=dynamics_ensemble,
            valid=valid,
            accepted=accepted,
            parent=rep_sample,
            details=sample_details,
            mover=self
        )

        move_details = MoveDetails()
        move_details.inputs = [rep_sample]
        move_details.trials = [trial]

#        new_set = SampleSet(samples=[sample], predecessor=globalstate, accepted=True)
#        new_set = globalstate.apply([sample], accepted = details.accepted, move=self)

        path = paths.SamplePathMoveChange(
                trials=[trial],
                mover=self,
                details=move_details
        )

        return path


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
            details.start_point.snapshot,
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
            details.start_point.snapshot.reversed,
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
        details.final_point = paths.ShootingPoint(self.selector, details.trial, len(partial_trajectory) - 1)

        pass


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
    def __init__(self, movers, ensembles=None,  weights=None):
        super(RandomChoiceMover, self).__init__(ensembles=ensembles)

        self.movers = movers

        if weights is None:
            self.weights = [1.0] * len(movers)
        else:
            self.weights = weights

        initialization_logging(init_log, self,
                               entries=['movers', 'weights'])

    @property
    def submovers(self):
        return self.movers

    @keep_selected_samples
    def move(self, globalstate):
        rand = np.random.random() * sum(self.weights)
        idx = 0
        prob = self.weights[0]
        while prob <= rand and idx < len(self.weights):
            idx += 1
            prob += self.weights[idx]

        logger_str = "{name} (RandomChoiceMover) selecting {mtype} (index {idx})"
        logger.info(logger_str.format(name=self.name, idx=idx, mtype=self.movers[idx].name))

        mover = self.movers[idx]

        details = MoveDetails()
        details.inputs = []
        details.choice = idx
        details.chosen_mover = mover

        path = paths.RandomChoicePathMoveChange(
            mover.move(globalstate),
            mover=self,
            details=details
        )

        return path


class ConditionalMover(PathMover):
    '''
    An if-then-else structure for PathMovers.

    Returns a SequentialPathMoveChange of the if_move movepath and the then_move
    movepath (if if_move is accepted) or the else_move movepath (if if_move
    is rejected).
    '''
    def __init__(self, if_mover, then_mover, else_mover, ensembles=None,
                 replicas='all'):
        super(ConditionalMover, self).__init__(ensembles=ensembles)
        self.if_mover = if_mover
        self.then_mover = then_mover
        self.else_mover = else_mover
        initialization_logging(init_log, self,
                               ['if_mover', 'then_mover', 'else_mover'])

    @property
    def submovers(self):
        return [self.if_mover, self.then_mover, self.else_mover]

    @keep_selected_samples
    def move(self, globalstate):
        subglobal = globalstate

        ifclause = self.if_mover.move(subglobal)
        samples = ifclause.samples
        subglobal = subglobal.apply_samples(samples)

        if ifclause.accepted:
            if self.then_mover is not None:
                resultclause = self.then_mover.move(subglobal)
            else:
                resultclause = paths.EmptyPathMoveChange()
        else:
            if self.else_mover is not None:
                resultclause = self.else_mover.move(subglobal)
            else:
                resultclause = paths.EmptyPathMoveChange()

        return paths.SequentialPathMoveChange([ifclause, resultclause], mover=self)


class SequentialMover(PathMover):
    '''
    Performs each of the moves in its movers list. Returns all samples
    generated, in the order of the mover list.

    For example, this would be used to create a move that does a sequence of
    replica exchanges in a given order, regardless of whether the moves
    succeed or fail.
    '''
    def __init__(self, movers, ensembles=None):
        super(SequentialMover, self).__init__(ensembles=ensembles)
        self.movers = movers
        initialization_logging(init_log, self, ['movers'])

    @property
    def submovers(self):
        return self.movers

    @keep_selected_samples
    def move(self, globalstate):
        logger.debug("Starting sequential move")

#        subglobal = SampleSet(self.legal_sample_set(globalstate))

        subglobal = globalstate
        pathmovechanges = []

        for mover in self.movers:
            logger.debug("Starting sequential move step "+str(mover))

            # Run the sub mover
            movepath = mover.move(subglobal)
            samples = movepath.samples
            subglobal = subglobal.apply_samples(samples)
            pathmovechanges.append(movepath)

        return paths.SequentialPathMoveChange(pathmovechanges, mover=self)


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
    @keep_selected_samples
    def move(self, globalstate):
        logger.debug("==== BEGINNING " + self.name + " ====")
        subglobal = paths.SampleSet(self.legal_sample_set(globalstate))
        pathmovechanges = []
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
            pathmovechanges.append(movepath)
            if not movepath.accepted:
                break

        logger.debug("==== FINISHING " + self.name + " ====")
        return paths.PartialAcceptanceSequentialPathMoveChange(pathmovechanges, mover=self)


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
    @keep_selected_samples
    def move(self, globalstate):
        logger.debug("Starting conditional sequential move")

#        subglobal = SampleSet(self.legal_sample_set(globalstate))

        subglobal = globalstate
        pathmovechanges = []

        for mover in self.movers:
            logger.debug("Starting sequential move step "+str(mover))

            # Run the sub mover
            movepath = mover.move(subglobal)
            samples = movepath.samples
            subglobal = subglobal.apply_samples(samples)
            pathmovechanges.append(movepath)

            if not movepath.accepted:
                break

        return paths.ConditionalSequentialPathMoveChange(pathmovechanges, mover=self)


class RestrictToLastSampleMover(PathMover):
    def __init__(self, mover):
        super(RestrictToLastSampleMover, self).__init__()
        self.mover = mover

    @keep_selected_samples
    def move(self, globalstate):
        movepath = self.mover.move(globalstate)
        return paths.KeepLastSamplePathMoveChange(movepath, mover=self)

    @property
    def submovers(self):
        return [self.mover]


class ReplicaIDChangeMover(PathMover):
    """
    Changes the replica ID for a path.
    """
    def __init__(self, replica_pairs, ensembles=None):
        self.replica_pairs = make_list_of_pairs(replica_pairs)
        super(ReplicaIDChangeMover, self).__init__(ensembles=ensembles)
        self._extra_details = ['rep_from', 'rep_to']
        initialization_logging(logger=init_log, obj=self,
                               entries=['replica_pairs'])

    @keep_selected_samples
    def move(self, globalstate):
        legal_from_rep = [rep[0] for rep in self.replica_pairs]
        rep_sample = self.select_sample(globalstate,
                                        ensembles=self.ensembles,
                                        replicas=legal_from_rep)

        legal_pairs = [pair for pair in self.replica_pairs
                       if pair[0]==rep_sample.replica]
        mypair = random.choice(legal_pairs)

        rep_from = mypair[0]
        rep_to = mypair[1]

        logger.info("Creating new sample from replica ID " + str(rep_from)
                    + " and putting it in replica ID " + str(rep_to))

        # note: currently this clones into a new replica ID. We might later
        # want to kill the old replica ID (and possibly rename this mover).

        sample_details = SampleDetails()

        new_sample = paths.Sample(
            replica=rep_to,
            ensemble=rep_sample.ensemble,
            trajectory=rep_sample.trajectory,
            valid=rep_sample.valid,
            accepted=True,
            parent=rep_sample,
            mover=self
        )

        # Can be used to remove the old sample. Not used yet!
        kill_sample = paths.Sample(
            replica=rep_from,
            trajectory=None,
            ensemble=rep_sample.ensemble,
            accepted=True,
            valid=True,
            parent=None,
            mover=self
        )

        details = MoveDetails()
        details.inputs = [rep_sample]
        details.trials = [rep_sample]
        details.mover = self
        setattr(details, 'rep_from', mypair[0])
        setattr(details, 'rep_to', mypair[1])

        return paths.SamplePathMoveChange(
            [new_sample],
            mover=self,
            details=details
        )


class EnsembleHopMover(PathMover):
    def __init__(self, bias=None, ensembles=None):
        # TODO: maybe allow a version of this with a single ensemble and ANY
        # ensemble can hop to that? messy to code; maybe same idea under
        # another name
        ensembles = make_list_of_pairs(ensembles)
        super(EnsembleHopMover, self).__init__(ensembles=ensembles)
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

    @keep_selected_samples
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

        valid = ens_to(trajectory)
        logger.info("Hop starts from legal ensemble: "+str(ens_from(trajectory)))
        logger.info("Hop ends in legal ensemble: "+str(ens_to(trajectory)))

        sample_details = SampleDetails()

        trial = paths.Sample(
            replica=replica,
            trajectory=trajectory,
            ensemble=ens_to,
            valid=valid,
            accepted=valid,
            details=sample_details,
            mover=self,
            parent=rep_sample
        )

        details = MoveDetails()
        details.inputs = [rep_sample]
        setattr(details, 'initial_ensemble', ens_from)
        setattr(details, 'trial_ensemble', ens_to)

        if valid:
            setattr(details, 'result_ensemble', ens_to)
        else:
            setattr(details, 'result_ensemble', ens_from)

        path = paths.SamplePathMoveChange(
            [trial],
            mover=self,
            details=details
        )

        return path


class ForceEnsembleChangeMover(EnsembleHopMover):
    '''
    Force an ensemble change in the sample.

    This should only be used as part of other moves, since this can create
    samples which are not valid.
    '''
    def __init__(self, ensembles=None):
        # no bias allowed
        super(ForceEnsembleChangeMover, self).__init__(ensembles=ensembles)

    @keep_selected_samples
    def move(self, globalstate):
        ens_pair = self.select_ensemble_pair(globalstate)
        ens_from = ens_pair[0]
        ens_to = ens_pair[1]
        rep_sample = self.select_sample(globalstate, ens_from)
        logger.debug("Selected sample: " + repr(rep_sample))

        replica = rep_sample.replica
        trajectory = rep_sample.trajectory

        sample_details = SampleDetails()

        sample = paths.Sample(
            trajectory=trajectory,
            ensemble=ens_to,
            replica=replica,
            details=sample_details,
            mover=self,
            parent=rep_sample
        )

        details = MoveDetails()
        details.accepted = True
        details.inputs = [rep_sample]
        setattr(details, 'initial_ensemble', ens_from)
        setattr(details, 'trial_ensemble', ens_to)
        setattr(details, 'result_ensemble', ens_to)

        path = paths.SamplePathMoveChange(
            [sample],
            mover=self,
            details=details
        )
        return path


class RandomSubtrajectorySelectMover(PathMover):
    '''
    Samples a random subtrajectory satisfying the given subensemble.

    If there are no subtrajectories which satisfy the subensemble, this
    returns the zero-length trajectory.
    '''
    def __init__(self, subensemble, n_l=None, ensembles=None):
        super(RandomSubtrajectorySelectMover, self).__init__(ensembles=ensembles)
        self.n_l=n_l
        self.subensemble = subensemble

    def _choose(self, trajectory_list):
        return random.choice(trajectory_list)

    @keep_selected_samples
    def move(self, globalstate):
        rep_sample = self.select_sample(globalstate)
        trajectory = rep_sample.trajectory
        replica = rep_sample.replica
        logger.debug("Working with replica " + str(replica) + " (" + str(trajectory) + ")")


        subtrajs = self.subensemble.split(trajectory)
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
            subtraj = paths.Trajectory([])

        sample_details = SampleDetails()

        sample = paths.Sample(
            replica=replica,
            trajectory=subtraj,
            ensemble=self.subensemble,
            valid=self.subensemble(subtraj),
            accepted=True,
            parent=rep_sample,
            mover=self
        )

        details = MoveDetails()
        details.inputs = [rep_sample]
        details.trials = [sample]

        path = paths.SamplePathMoveChange(
            [sample],
            mover=self,
            details=details
        )

        return path


class FirstSubtrajectorySelectMover(RandomSubtrajectorySelectMover):
    '''
    Samples the first subtrajectory satifying the given subensemble.

    If there are no subtrajectories which satisfy the ensemble, this returns
    the zero-length trajectory.
    '''
    def _choose(self, trajectory_list):
        return trajectory_list[0]


class FinalSubtrajectorySelectMover(RandomSubtrajectorySelectMover):
    '''
    Samples the final subtrajectory satifying the given subensemble.

    If there are no subtrajectories which satisfy the ensemble, this returns
    the zero-length trajectory.
    '''
    def _choose(self, trajectory_list):
        return trajectory_list[-1]


class PathReversalMover(PathMover):

    @keep_selected_samples
    def move(self, globalstate):
        rep_sample = self.select_sample(globalstate, self.ensembles)

        trajectory = rep_sample.trajectory
        ensemble = rep_sample.ensemble
        replica = rep_sample.replica


        reversed_trajectory = trajectory.reversed

        valid = ensemble(reversed_trajectory)
        logger.info("PathReversal move accepted: "+str(valid))

        acceptance_probability = 0.0

        if valid:
            acceptance_probability = 1.0

        sample_details = SampleDetails()
        sample_details.acceptance_probability = acceptance_probability

        trial = paths.Sample(
            replica=replica,
            trajectory=reversed_trajectory,
            ensemble=ensemble,
            valid=valid,
            accepted=valid,
            details=sample_details,
            mover=self,
            parent=rep_sample
        )

        details = MoveDetails()
        details.inputs = [rep_sample]
        details.trials = [trial]

        path = paths.SamplePathMoveChange(
            [trial],
            mover=self,
            details=details
        )

        return path


class ReplicaExchangeMover(PathMover):
    def __init__(self, bias=None, ensembles=None):
        replicas = 'all'
        ensembles = make_list_of_pairs(ensembles)
        # either replicas or ensembles must be a list of pairs; more
        # complicated filtering can be done with a wrapper class
        super(ReplicaExchangeMover, self).__init__(ensembles=ensembles)
        # TODO: add support for bias; cf EnsembleHopMover
        self.bias = bias
        initialization_logging(logger=init_log, obj=self,
                               entries=['bias'])

    @keep_selected_samples
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
        logger.info("trajectory " + repr(trajectory1) +
                    " into ensemble " + repr(ensemble2) +
                    " : " + str(from1to2))
        from2to1 = ensemble1(trajectory2)
        logger.info("trajectory " + repr(trajectory2) +
                    " into ensemble " + repr(ensemble1) +
                    " : " + str(from2to1))

        accepted = from1to2 and from2to1

        trial1 = paths.Sample(
            replica=replica1,
            trajectory=trajectory1,
            ensemble=ensemble2,
            valid=from1to2,
            accepted=accepted,
            parent=s1,
            details = SampleDetails(),
            mover=self
        )
        trial2 = paths.Sample(
            replica=replica2,
            trajectory=trajectory2,
            ensemble=ensemble1,
            valid=from2to1,
            accepted=accepted,
            parent=s2,
            details=SampleDetails(),
            mover=self
        )

        details = MoveDetails()
        details.inputs = [s1, s2]
        details.trials = [trial1, trial2]
        setattr(details, 'ensembles', [ensemble1, ensemble2])

        path = paths.SamplePathMoveChange(
            [trial1, trial2],
            mover=self,
            details=details
        )

        return path


# TODO: Turn Filter into real mover with own movechange ?

class FilterByReplica(PathMover):
    def __init__(self, mover, replicas):
        if type(replicas) is not list:
            replicas = [replicas]
        self.replicas = replicas
        self.mover = mover
        # TODO: clean this up
        pass

    @keep_selected_samples
    def move(self, globalstate):
        filtered_gs = paths.SampleSet(
            [s for s in globalstate if s.replica in self.replicas]
        )
        return self.mover.move(filtered_gs)


class FilterBySample(PathMover):
    def __init__(self, mover, selected_samples, use_all_samples=None):
        if type(selected_samples) is not list:
            selected_samples = [selected_samples]
        self.selected_samples = selected_samples
        self.mover = mover
        self.use_all_samples = use_all_samples

    @keep_selected_samples
    def move(self, globalstate):
        return paths.FilterSamplesPathMoveChange(
            self.mover.move(globalstate),
            selected_samples=self.selected_samples,
            use_all_samples=self.use_all_samples,
            mover=self
        )


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
    def __init__(self, selector, ensembles=None):
        movers = [
            ForwardShootMover(selector, ensembles),
            BackwardShootMover(selector, ensembles)
        ]
        super(OneWayShootingMover, self).__init__(
            movers=movers, ensembles=ensembles
        )
        self.selector = selector

class MinusMover(ConditionalSequentialMover, ReplicaExchangeMover):
    '''
    Instance of a MinusMover.

    The minus move combines a replica exchange with path extension to swap
    paths between the innermost regular TIS interface ensemble and the minus
    interface ensemble. This is particularly useful for improving sampling
    of path space.

    Note that the inheritance from ReplicaExchangeMover is only to assist
    with `isinstance` in later analysis. Since the only two functions here
    are `.__init__() and `.move()`, both of which exist in both parent
    classes, the calls to `super` will use the version in
    ConditionalSequentalMover. However, analysis routines will see
    `isinstance(minus, ReplicaExchangeMover)` as True.
    '''
    def __init__(self, minus_ensemble, innermost_ensemble, 
                 ensembles=None):
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

        super(MinusMover, self).__init__(movers=movers, ensembles=ensembles)


    @keep_selected_samples
    def move(self, globalstate):
        result = super(MinusMover, self).move(globalstate).closed
        return result


class PathSimulatorMover(PathMover):
    """
    This just wraps a mover and references the used pathsimulator
    """
    def __init__(self, mover, pathsimulator):
        super(PathSimulatorMover, self).__init__()
        self.mover = mover
        self.pathsimulator = pathsimulator

    def move(self, globalstate, step=-1):
        return paths.PathSimulatorPathMoveChange(
            self.mover.move(globalstate),
            self.pathsimulator,
            step=step,
            mover=self
        )

    @property
    def submovers(self):
        return [self.mover]


class MultipleSetMinusMover(RandomChoiceMover):
    pass

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

        mover_set = []
        for (selector, iface) in zip(selector_set, interface_set):
            mover = OneWayShootingMover(selector=selector,
                                        ensembles=[iface])
            mover.name = "OneWayShootingMover " + str(iface.name)
            mover_set.append(mover)

        return mover_set

    @staticmethod
    def TwoWayShootingSet():
        pass

    @staticmethod
    def NearestNeighborRepExSet():
        pass



class Details(OPSObject):
    '''Details of an object. Can contain any data
    '''

    def __init__(self, **kwargs):
        for key, value in kwargs:
            setattr(self, key, value)

    def __str__(self):
        # primarily for debugging/interactive use
        mystr = ""
        for key in self.__dict__.keys():
            if not isinstance(self.__dict__[key], paths.Ensemble):
                mystr += str(key) + " = " + str(self.__dict__[key]) + '\n'
        return mystr


class MoveDetails(Details):
    '''Details of the move as applied to a given replica

    Attributes
    ----------
    replica : integer
        replica ID to which this trial move would apply
    inputs : list of Trajectory
        the Samples which were used as inputs to the move
    trial : Trajectory
        the Trajectory
    trial_is_in_ensemble : bool
        whether the attempted move created a trajectory in the right
        ensemble
    mover : PathMover
        the PathMover which generated this sample out of other samples

    Specific move types may have add several other attributes for each
    MoveDetails object. For example, shooting moves will also include
    information about the shooting point selection, etc.

    TODO (or at least to put somewhere):
    rejection_reason : String
        explanation of reasons the path was rejected

    RENAME: inputs=>initial
            accepted=>trial_in_ensemble (probably only in shooting)

    TODO:
    Currently trial/accepted are in terms of Trajectory objects. I
    think it makes more sense for them to be Samples.
    I kept trial, accepted as a trajectory and only changed inputs
    to a list of samples. Since trial, accepted are move related
    to the shooting and not necessarily dependent on a replica or
    initial ensemble.
    '''

    def __init__(self, **kwargs):
        self.inputs=None
        self.trials=None
        self.results=None
        super(MoveDetails, self).__init__(**kwargs)

    # @staticmethod
    # def initialization(sample):
    #     return MoveDetails.initialization_from_scratch(sample.trajectory,
    #                                                    sample.ensemble)
    #
    # @staticmethod
    # def initialization_from_scratch(trajectory, ensemble):
    #     details = MoveDetails()
    #     details.inputs = []
    #     details.trials = trajectory
    #     details.ensemble = ensemble # might go to Change and not Details
    #     details.results = trajectory
    #     return details


class SampleDetails(Details):
    '''Details of a sample

    Attributes
    ----------
    selection_probability : float
        the chance that a sample will be accepted due to asymmetrical proposal
    '''

    def __init__(self, **kwargs):
        self.selection_probability=1.0
        super(SampleDetails, self).__init__(**kwargs)
