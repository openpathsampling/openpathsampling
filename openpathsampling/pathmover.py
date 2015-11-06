"""
Created on 19.07.2014

@author: Jan-Hendrik Prinz
@author: David W. H. Swenson
"""

import random
import logging

import numpy as np

import openpathsampling as paths
from openpathsampling.base import StorableNamedObject, StorableObject
from ops_logging import initialization_logging
from treelogic import TreeMixin


logger = logging.getLogger(__name__)
init_log = logging.getLogger('openpathsampling.initialization')

# TODO: Remove if really not used anymore otherwise might move to utils or tools
def make_list_of_pairs(l):
    """
    Converts input from several possible formats into a list of pairs: used
    to clean input for swap-like moves.

    Allowed input formats:
    * flat list of length 2N
    * list of pairs
    * None (returns None)

    Anything else will lead to a ValueError or AssertionError
    """
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

    if list_of_lists:
        for elem in l:
            assert len(elem)==2, "List of lists: inner list length != 2"
        outlist = l
    else:
        assert len(l) % 2 == 0, "Flattened list: length not divisible by 2"
        outlist = [
            [a, b] for (a, b) in zip(l[slice(0, None, 2)], l[slice(1, None,2 )])
        ]
    # Note that one thing we don't check is whether the items are of the
    # same type. That might be worth doing someday; for now, we trust that
    # part to work.
    return outlist


class PathMover(TreeMixin, StorableNamedObject):
    """
    A PathMover is the description of a move in replica space.
    
    Notes
    -----
    A pathmover takes a SampleSet() and returns PathMoveChange() that is
    used to change the old SampleSet() to the new one.

    SampleSet1 + PathMoveChange1 => SampleSet2

    A PathMoveChange is effectively a list of Samples. The change acts upon
    a SampleSet by replacing existing Samples in the same ensemble
    sequentially.

    SampleSet({samp1(ens1), samp2(ens2), samp3(ens3)}) +
        PathMoveChange([samp4(ens2)])
        => SampleSet({samp1(ens1), samp4(ens2), samp3(ens3)})

    Note, that a SampleSet is an unordered list (or a set). Hence the ordering
    in the example is arbitrary.
    
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
    """

    def __init__(self):
        StorableNamedObject.__init__(self)

        self._in_ensembles = None
        self._out_ensembles = None
        self._len = None

#        initialization_logging(logger=init_log, obj=self,
#                               entries=['ensembles'])

    _is_ensemble_change_mover = None
    @property
    def is_ensemble_change_mover(self):
        if self._is_ensemble_change_mover is None:
            return False
        else:
            return self._is_ensemble_change_mover

    _is_canonical = None
    @property
    def is_canonical(self):
        return self._is_canonical


    @property
    def default_name(self):
        return self.__class__.__name__[:-5]

    # +-------------------------------------------------------------------------
    # | tree implementation overrides
    # +-------------------------------------------------------------------------

    @property
    def _subnodes(self):
        return self.submovers

    @property
    def identifier(self):
        return self

    @staticmethod
    def _default_match(original, test):
        if isinstance(test, paths.PathMover):
            return original is test
        elif issubclass(test, paths.PathMover):
            return original.__class__ is test
        else:
            return False

    @property
    def submovers(self):
        """
        Returns a list of submovers

        Returns
        -------
        list of openpathsampling.PathMover
            the list of sub-movers
        """
        return []

    @staticmethod
    def _flatten(ensembles):
        if type(ensembles) is list:
            return [s for ens in ensembles for s in PathMover._flatten(ens)]
        else:
            return [ensembles]

    def _ensemble_signature(self, as_set=False):
        """Return tuple form of (input_ensembles, output_ensembles).
        
        Useful for MoveScheme, e.g., identifying which movers should be
        removed as part of a replacement.
        """
        inp = tuple(self.input_ensembles)
        out = tuple(self.output_ensembles)
        if as_set:
            inp = set(inp)
            out = set(out)
        return (inp, out)
               
    @property
    def ensemble_signature(self):
        return self._ensemble_signature(as_set=False)

    @property
    def ensemble_signature_set(self):
        return self._ensemble_signature(as_set=True)

    @property
    def input_ensembles(self):
        """Return a list of possible used ensembles for this mover

        This list contains all Ensembles from which this mover might pick
        samples. This is very useful to determine on which ensembles a
        mover acts for analysis and sanity checking.

        Returns
        -------
        list of Ensemble
            the list of input ensembles
        """
        if self._in_ensembles is None:
            ensembles = self._get_in_ensembles()

            self._in_ensembles = list(set(self._flatten(ensembles)))

        return self._in_ensembles

    @property
    def output_ensembles(self):
        """Return a list of possible returned ensembles for this mover

        This list contains all Ensembles for which this mover might return
        samples. This is very useful to determine on which ensembles a
        mover affects in later steps for analysis and sanity checking.

        Returns
        -------
        list of Ensemble
            the list of output ensembles
        """

        if self._out_ensembles is None:
            ensembles = self._get_out_ensembles()

            self._out_ensembles = list(set(self._flatten(ensembles)))

        return self._out_ensembles

    def _get_in_ensembles(self):
        """Function that computes the list of input ensembles
        """
        return []

    def _get_out_ensembles(self):
        """Function that computes the list of output ensembles

        Default is the same as in_ensembles
        """
        return self._get_in_ensembles()

    def legal_sample_set(self, globalstate, ensembles=None, replicas='all'):
        """
        This returns all the samples from globalstate which are in both
        self.replicas and the parameter ensembles. If ensembles is None, we
        use self.ensembles. If you want all ensembles allowed, pass
        ensembles='all'.
        """
        mover_replicas = globalstate.replica_list()

        if replicas == 'all':
            selected_replicas = globalstate.replica_list()
        else:
            selected_replicas = replicas

        reps = list(set(mover_replicas) & set(selected_replicas))
        rep_samples = []
        for rep in reps:
            rep_samples.extend(globalstate.all_from_replica(rep))

        # logger.debug("ensembles = " + str([ensembles]))
        # logger.debug("self.ensembles = " + str(self.ensembles))
        if ensembles is None:
            ensembles = 'all'

        if ensembles == 'all':
            legal_samples = rep_samples
        else:
            ens_samples = []
            if type(ensembles) is not list:
                ensembles = [ensembles]
            for ens in ensembles:
                # try:
                #     ens_samples.extend(globalstate.all_from_ensemble(ens[0]))
                # except TypeError:
                ens_samples.extend(globalstate.all_from_ensemble(ens))
            legal_samples = list(set(rep_samples) & set(ens_samples))

        return legal_samples

    def select_sample(self, globalstate, ensembles=None, replicas=None):
        """
        Returns one of the legal samples given self.replica and the ensemble
        set in ensembles.

        TODO: This must be saved somehow (it is actually I think), otherwise
        Samples are not reproducible when applied to a SampleSet!
        """
        if replicas is None:
            replicas = 'all'

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
        """
        Run the generation starting with the initial globalstate specified.

        Parameters
        ----------
        globalstate : SampleSet
            the initially used sampleset
        
        Returns
        -------        
        samples : PathMoveChange
            the PathMoveChange instance describing the change from the old to
            the new SampleSet

        """

        return paths.EmptyPathMoveChange()  # pragma: no cover

    def __str__(self):
        if self.name == self.__class__.__name__:
            return self.__repr__()
        else:
            return self.name


###############################################################################
# MOVER TYPES
###############################################################################

class MoverType(object):
    pass


# TODO: empty class. Remove
class SwappingMover(MoverType):
    """
    A mover that swaps samples from ensembles in some way. Relevant for mixing
    """


###############################################################################
# GENERATORS
###############################################################################

class SampleMover(PathMover):
    engine = None

    def __init__(self):
        super(SampleMover, self).__init__()

    @classmethod
    def metropolis(cls, trials):
        """Implements the Metropolis acceptance for a list of trial samples

        The Metropolis uses the .bias for each sample and checks of samples
        are valid - are in the proposed ensemble. This will give an acceptance
        probability for all samples. If the product is smaller than a random
        number the change will be accepted.

        Parameters
        ----------
        trials : list of openpathsampling.Sample
            the list of all samples to be applied in a change.

        Returns
        -------
        bool
            True if the trial is accepted, False otherwise
        details : openpathsampling.MoveDetails
            Returns a MoveDetails object that contains information about the
            decision, i.e. total acceptance and random number

        """

        shoot_str = "MC in {cls} using samples {trials}"
        logger.info(shoot_str.format(cls=cls.__name__, trials=trials))
        trial_dict = dict()
        for trial in trials:
            trial_dict[trial.ensemble] = trial

        accepted = True
        probability = 1.0

        # TODO: This isn't right. `bias` should be associated with the 
        # change; not with each individual sample. ~~~DWHS
        for ens, sample in trial_dict.iteritems():
            valid = ens(sample.trajectory)
            if not valid:
                # one sample not valid reject
                accepted = False
                probability = 0.0
                break
            else:
                probability *= sample.bias

        rand = random.random()

        if rand > probability:
            # rejected
            accepted = False

        details = paths.MoveDetails(
            total_acceptance=probability,
            random_value=rand
        )

        return accepted, details

    @property
    def submovers(self):
        # Movers do not have submovers!
        return []

    def _called_ensembles(self):
        """Function to determine which ensembles to pick samples from

        Returns
        -------
        list of Ensemble
            the list of ensembles. Samples can then be selected using
            PathMover.select_sample
        """

        # Default is that the list of ensembles is in self.ensembles
        return []

    def move(self, globalstate):
        # 1. pick a set of ensembles (in case we allow to pick several ones)
        ensembles = self._called_ensembles()

        # 2. pick samples from these ensembles
        samples = [self.select_sample(globalstate, ens) for ens in ensembles]

        # 3. pass these samples to the generator
        trials = self(*samples)

        # 4. accept/reject
        accepted, details = self._accept(trials)

        # 5. and return a PMC
        if accepted:
            return paths.AcceptedSamplePathMoveChange(
                samples=trials,
                mover=self,
                details=details
            )
        else:
            return paths.RejectedSamplePathMoveChange(
                samples=trials,
                mover=self,
                details=details
            )

    def __call__(self, *args):
        """Generate trial samples directly

        PathMovers can also be called directly with a list of samples that are
        then used to generate new samples. If the Mover is used as a move
        the move will first determine the input samples and then pass these to
        this function
        """

        # Default is that the original samples are returned
        return args

    def _accept(self, trials):
        """Function to determine the acceptance of a trial

        Defaults to calling the Metropolis acceptance criterion for all returned
        trial samples. Means all samples most be valid and accepted.
        """
        return self.metropolis(trials)


###############################################################################
# SHOOTING GENERATORS
###############################################################################

class EngineMover(SampleMover):
    """Baseclass for Movers that use an engine
    """
    engine = None
    def __init__(self, ensemble, target_ensemble, selector):
        super(EngineMover, self).__init__()
        self.selector = selector
        self.ensemble = ensemble
        self.target_ensemble = target_ensemble

    def _called_ensembles(self):
        return [self.ensemble]

    def _get_in_ensembles(self):
        return [self.ensemble]

    def _get_out_ensembles(self):
        return [self.target_ensemble]

    def __call__(self, input_sample):
        initial_trajectory = input_sample.trajectory
        replica = input_sample.replica

        shooting_index = self.selector.pick(initial_trajectory)

        trial_trajectory = self._run(initial_trajectory, shooting_index)


        bias = self.selector.probability_ratio(
            initial_trajectory[shooting_index],
            initial_trajectory,
            trial_trajectory
        )

        # temporary test to make sure nothing went weird
        # old_bias = initial_point.sum_bias / trial_point.sum_bias
        # assert(abs(bias - old_bias) < 10e-6)
        assert(initial_trajectory[shooting_index] in trial_trajectory)

        # we need to save the initial
        trial_details = paths.SampleDetails(
            initial_trajectory=initial_trajectory,
            shooting_snapshot=initial_trajectory[shooting_index]
        )

        trial = paths.Sample(
            replica=replica,
            trajectory=trial_trajectory,
            ensemble=self.target_ensemble,
            parent=input_sample,
            details=trial_details,
            mover=self,
            bias=bias
        )

        trials = [trial]

        return trials

    def _make_forward_trajectory(self, trajectory, shooting_index):
        initial_snapshot = trajectory[shooting_index].copy()
        run_f = paths.PrefixTrajectoryEnsemble(self.target_ensemble, 
                                               trajectory[0:shooting_index]
                                              ).can_append
        partial_trajectory = self.engine.generate(initial_snapshot, 
                                                  running=[run_f])
        # keep the original snapshot in the trial_trajectory
        trial_trajectory = (trajectory[0:shooting_index + 1] 
                            + partial_trajectory[1:])
        return trial_trajectory

    def _make_backward_trajectory(self, trajectory, shooting_index):
        initial_snapshot = trajectory[shooting_index].reversed_copy()
        run_f = paths.SuffixTrajectoryEnsemble(self.target_ensemble,
                                               trajectory[shooting_index + 1:]
                                              ).can_prepend
        partial_trajectory = self.engine.generate(initial_snapshot, 
                                                  running=[run_f])
        # keep the original snapshot in the trial_trajectory
        trial_trajectory = (partial_trajectory.reversed[:-1] +
                            trajectory[shooting_index:])
        return trial_trajectory


    def _run(self, trajectory, shooting_index):
        """Takes initial trajectory and shooting point; return trial
        trajectory"""
        shoot_str = "Running {sh_dir} from frame {fnum} in [0:{maxt}]"
        logger.info(shoot_str.format(
            fnum=shooting_index,
            maxt=len(trajectory)-1,
            sh_dir=self._direction
        ))

        if self._direction == "forward":
            trial_trajectory = self._make_forward_trajectory(
                trajectory, shooting_index
            )
        elif self._direction == "backward":
            trial_trajectory = self._make_backward_trajectory(
                trajectory, shooting_index
            )
        else:
            raise RuntimeError("Unknown direction: " + str(self._direction))

        return trial_trajectory


class ForwardShootMover(EngineMover):
    """A forward shooting sample generator
    """
    _direction = "forward"
    def __init__(self, ensemble, selector):
        super(ForwardShootMover, self).__init__(
            ensemble=ensemble,
            target_ensemble=ensemble,
            selector=selector
        )


class BackwardShootMover(EngineMover):
    """A Backward shooting generator
    """
    _direction = "backward"
    def __init__(self, ensemble, selector):
        super(BackwardShootMover, self).__init__(
            ensemble=ensemble,
            target_ensemble=ensemble,
            selector=selector
        )


class ForwardExtendMover(EngineMover):
    """
    A Sample Mover implementing Forward Extension
    """
    _direction = "forward"
    def __init__(self, ensemble, target_ensemble):
        super(ForwardExtendMover, self).__init__(
            ensemble=ensemble,
            target_ensemble=target_ensemble,
            selector=paths.FinalFrameSelector(),
        )


class BackwardExtendMover(EngineMover):
    """
    A Sample Mover implementing Backward Extension
    """
    _direction = "backward"
    def __init__(self, ensemble, target_ensemble):
        super(BackwardExtendMover, self).__init__(
            ensemble=ensemble,
            target_ensemble=target_ensemble,
            selector=paths.FirstFrameSelector(),
        )


###############################################################################
# REPLICA EXCHANGE GENERATORS
###############################################################################

class ReplicaExchangeMover(SampleMover):
    """
    A Sample Mover implementing a standard Replica Exchange
    """
    _is_ensemble_change_mover = True

    def __init__(self, ensemble1, ensemble2, bias=None):
        """
        Parameters
        ----------
        ensemble1 : openpathsampling.Ensemble
            one of the ensemble between to make the repex move
        ensemble2 : openpathsampling.Ensemble
            one of the ensemble between to make the repex move
        bias : list of float
            bias is not used yet

        """
        # either replicas or ensembles must be a list of pairs; more
        # complicated filtering can be done with a wrapper class
        super(ReplicaExchangeMover, self).__init__()
        # TODO: add support for bias; cf EnsembleHopMover
        self.bias = bias
        self.ensemble1 = ensemble1
        self.ensemble2 = ensemble2

        initialization_logging(logger=init_log, obj=self,
                               entries=['bias', 'ensemble1', 'ensemble2'])

    def _called_ensembles(self):
        return [self.ensemble1, self.ensemble2]

    def _get_in_ensembles(self):
        return [self.ensemble1, self.ensemble2]

    def _get_out_ensembles(self):
        return [self.ensemble1, self.ensemble2]

    def __call__(self, sample1, sample2):
        # convert sample to the language used here before
        trajectory1 = sample1.trajectory
        trajectory2 = sample2.trajectory
        ensemble1 = sample1.ensemble
        ensemble2 = sample2.ensemble
        replica1 = sample1.replica
        replica2 = sample2.replica

        from1to2 = ensemble2(trajectory1)
        logger.debug("trajectory " + repr(trajectory1) +
                     " into ensemble " + repr(ensemble2) +
                     " : " + str(from1to2))
        from2to1 = ensemble1(trajectory2)
        logger.debug("trajectory " + repr(trajectory2) +
                     " into ensemble " + repr(ensemble1) +
                     " : " + str(from2to1))

        trial1 = paths.Sample(
            replica=replica1,
            trajectory=trajectory1,
            ensemble=ensemble2,
            parent=sample1,
            details=SampleDetails(),
            mover=self
        )
        trial2 = paths.Sample(
            replica=replica2,
            trajectory=trajectory2,
            ensemble=ensemble1,
            parent=sample2,
            details=SampleDetails(),
            mover=self
        )

        return [trial1, trial2]


class StateSwapMover(SampleMover):
    def __init__(self, ensemble1, ensemble2, bias=None):
        """
        A move to swap states for state changing smaples

        This does a replica exchange with prededing PathReversal and
        will only succeed if initial and final state are different

        Parameters
        ----------
        ensemble1 : openpathsampling.Ensemble
            one of the ensemble between to make the swap move
        ensemble2 : openpathsampling.Ensemble
            one of the ensemble between to make the swap move
        bias : list of float
            bias is not used yet

        Notes
        -----
        So, if ensemble1 goes from A to B, then ensemble2 must go from B to A.
        """
        # either replicas or ensembles must be a list of pairs; more
        # complicated filtering can be done with a wrapper class
        super(StateSwapMover, self).__init__()
        self.bias = bias
        self.ensemble1 = ensemble1
        self.ensemble2 = ensemble2

        initialization_logging(logger=init_log, obj=self,
                               entries=['bias', 'ensemble1', 'ensemble2'])

    def _called_ensembles(self):
        return [self.ensemble1, self.ensemble2]

    def _get_in_ensembles(self):
        return [self.ensemble1, self.ensemble2]

    def _get_out_ensembles(self):
        return [self.ensemble1, self.ensemble2]

    def __call__(self, sample1, sample2):
        # convert sample to the language used here before

        # it is almost a RepEx move but the two trajectories are reversed
        trajectory1 = sample1.trajectory.reversed
        trajectory2 = sample2.trajectory.reversed
        ensemble1 = sample1.ensemble
        ensemble2 = sample2.ensemble
        replica1 = sample1.replica
        replica2 = sample2.replica

        from1to2 = ensemble2(trajectory1)
        logger.debug("trajectory " + repr(trajectory1) +
                     " into ensemble " + repr(ensemble2) +
                     " : " + str(from1to2))
        from2to1 = ensemble1(trajectory2)
        logger.debug("trajectory " + repr(trajectory2) +
                     " into ensemble " + repr(ensemble1) +
                     " : " + str(from2to1))

        trial1 = paths.Sample(
            replica=replica1,
            trajectory=trajectory1,
            ensemble=ensemble2,
            parent=sample1,
            details=SampleDetails(),
            mover=self
        )
        trial2 = paths.Sample(
            replica=replica2,
            trajectory=trajectory2,
            ensemble=ensemble1,
            parent=sample2,
            details=SampleDetails(),
            mover=self
        )

        return [trial1, trial2]


###############################################################################
# SUBTRAJECTORY GENERATORS
###############################################################################


class RandomSubtrajectorySelectMover(SampleMover):
    """
    Samples a random subtrajectory satisfying the given subensemble.

    If there are no subtrajectories which satisfy the subensemble, this
    returns the zero-length trajectory.

    Parameters
    ----------
    ensemble : openpathsampling.Ensemble
        the set of allows samples to chose from
    subensemble : openpathsampling.Ensemble
        the subensemble to be searched for
    n_l : int or None
        the number of subtrajectories that need to be found. If
        `None` every number of subtrajectories > 0 is okay.
        Otherwise the move is only accepted if exactly n_l subtrajectories
        are found.

    """
    _is_ensemble_change_mover = True
    def __init__(self, ensemble, sub_ensemble, n_l=None):
        super(RandomSubtrajectorySelectMover, self).__init__(
        )
        self.n_l = n_l
        self.ensemble = ensemble
        self.sub_ensemble = sub_ensemble

    def _called_ensembles(self):
        return [ self.ensemble ]

    def _get_in_ensembles(self):
        return [ self.ensemble ]

    def _get_out_ensembles(self):
        return [ self.sub_ensemble ]

    def _choose(self, trajectory_list):
        return random.choice(trajectory_list)

    def __call__(self, trial):
        initial_trajectory = trial.trajectory
        replica = trial.replica
        logger.debug("Working with replica " + str(replica) + " (" + str(initial_trajectory) + ")")

        subtrajs = self.sub_ensemble.split(initial_trajectory)
        logger.debug("Found "+str(len(subtrajs))+" subtrajectories.")

        if (self.n_l is None and len(subtrajs) > 0) or \
            (self.n_l is not None and len(subtrajs) == self.n_l):
            subtraj = self._choose(subtrajs)

            bias = 1.0

            trial = paths.Sample(
                replica=replica,
                trajectory=subtraj,
                ensemble=self.sub_ensemble,
                parent=trial,
                mover=self,
                bias=bias
            )

            trials = [trial]
        else:
            trials = []

        return trials



class FirstSubtrajectorySelectMover(RandomSubtrajectorySelectMover):
    """
    Samples the first subtrajectory satifying the given subensemble.

    If there are no subtrajectories which satisfy the ensemble, this returns
    the zero-length trajectory.
    """
    def _choose(self, trajectory_list):
        return trajectory_list[0]


class FinalSubtrajectorySelectMover(RandomSubtrajectorySelectMover):
    """
    Samples the final subtrajectory satifying the given subensemble.

    If there are no subtrajectories which satisfy the ensemble, this returns
    the zero-length trajectory.
    """
    def _choose(self, trajectory_list):
        return trajectory_list[-1]

###############################################################################
# REVERSAL GENERATOR
###############################################################################

class PathReversalMover(SampleMover):

    def __init__(self, ensemble):
        """
        Parameters
        ----------
        ensemble : openpathsampling.Ensemble
            the specific ensemble to be reversed in
        """
        super(PathReversalMover, self).__init__()
        self.ensemble = ensemble

    def _called_ensembles(self):
        return [ self.ensemble ]

    def _get_in_ensembles(self):
        return [ self.ensemble ]

    def __call__(self, trial):
        trajectory = trial.trajectory
        ensemble = trial.ensemble
        replica = trial.replica

        reversed_trajectory = trajectory.reversed

        valid = ensemble(reversed_trajectory)
        logger.info("PathReversal move accepted: "+str(valid))

        bias = 1.0

        trial = paths.Sample(
            replica=replica,
            trajectory=reversed_trajectory,
            ensemble=ensemble,
            mover=self,
            parent=trial,
            bias=bias
        )

        return [trial]


class EnsembleHopMover(SampleMover):
    _is_ensemble_change_mover = True
    def __init__(self, ensemble, target_ensemble, change_replica=None, bias=None):
        """
        Parameters
        ----------
        ensemble : openpathsampling.Ensemble
            the initial ensemble to be jumped from
        target_ensemble : openpathsampling.Ensemble
            the final ensemble to be jumped to
        change_replica : int of None
            if None the replica id of the chosen sample will not be changed.
            Otherwise the replica id will be set to change_replica. This is
            useful when hoping to ensembles to create a new replica.
        bias : float, dict or None (default)
            gives the bias of accepting (not proposing) a hop. A float will
            be the acceptance for all possible attempts. If a dict is given,
            then it contains a list of ensembles and a matrix. None means
            no bias

        Notes
        -----
        The bias dict has the following form :
            { 'ensembles' : [ens_1, ens_2, ens_n],
              'values' : np.array((n,n))
              }
        The numpy array contains all the acceptance probabilties. If possible
        a HopMover should (as all movers) be used for only a specific hop and
        not multiple ones.
        """
        super(EnsembleHopMover, self).__init__()
        # ensembles -- another version might take a value for each ensemble,
        # and use the ratio; this latter is better for CITIS
        self.ensemble = ensemble
        self.target_ensemble = target_ensemble
        self.bias = bias
        self.change_replica = change_replica

        initialization_logging(
            logger=init_log,
            obj=self,
            entries=['bias']
        )

    def _called_ensembles(self):
        return [ self.ensemble ]

    @property
    def submovers(self):
        return []

    def _get_in_ensembles(self):
        return [self.ensemble]

    def _get_out_ensembles(self):
        return [self.target_ensemble]

    def __call__(self, rep_sample):
        ens_from = self.ensemble
        ens_to = self.target_ensemble

        logger.debug("Selected sample: " + repr(rep_sample))
        replica = rep_sample.replica

        if self.change_replica is not None:
            replica = self.change_replica

        logger.info("Attempting ensemble hop from {e1} to {e2} replica ID {rid}".format(
            e1=repr(ens_from), e2=repr(ens_to), rid=repr(replica)))

        trajectory = rep_sample.trajectory
        logger.debug("  selected replica: " + str(replica))
        logger.debug("  initial ensemble: " + repr(rep_sample.ensemble))

        logger.info("Hop starts from legal ensemble: "+str(ens_from(trajectory)))
        logger.info("Hop ends in legal ensemble: "+str(ens_to(trajectory)))

        sample_details = SampleDetails()

        if type(self.bias) is float:
            bias = self.bias
        elif type(self.bias) is dict:
            # special dict
            ens = self.bias['ensembles']
            e1 = ens.index(ens_from)
            e2 = ens.index(ens_to)
            bias = float(self.bias['values'][e1,e2])
        else:
            bias = 1.0

        trial = paths.Sample(
            replica=replica,
            trajectory=trajectory,
            ensemble=ens_to,
            details=sample_details,
            mover=self,
            parent=rep_sample,
            bias=bias
        )

        details = MoveDetails()
        setattr(details, 'initial_ensemble', ens_from)
        setattr(details, 'trial_ensemble', ens_to)
        setattr(details, 'bias', bias)

        return [trial]


# ****************************************************************************
#  SELECTION MOVERS
# ****************************************************************************

class SelectionMover(PathMover):
    """
    A general mover that selects a single mover from a set of possibilities

    This is a basic class for all sorts of selectors, like RandomChoice,
    RandomAllowedChoice. The way it works is to generate a list of weights
    and pick a random one using the weights. This is as general as possible
    and is chosen because it also allows to store the possibilities in a
    general way for better comparison

    Attributes
    ----------
    movers : list of openpathsampling.PathMover
        the PathMovers to choose from
    """

    def __init__(self, movers):
        super(SelectionMover, self).__init__()

        self.movers = movers

        initialization_logging(init_log, self,
                               entries=['movers'])

    @property
    def submovers(self):
        return self.movers

    @property
    def is_ensemble_change_mover(self):
        if self._is_ensemble_change_mover is not None:
            return self._is_ensemble_change_mover
        sub_change = False
        for mover in self.movers:
            if mover.is_ensemble_change_mover:
                sub_change = True
                break
        return sub_change


    def _get_in_ensembles(self):
        return [ sub.input_ensembles for sub in self.submovers ]

    def _get_out_ensembles(self):
        return [ sub.output_ensembles for sub in self.submovers ]

    def _selector(self, globalstate):
        # Default always picks by random choice
        return [1.0] * len(self.movers)

    def move(self, globalstate):
        weights = self._selector(globalstate)

        rand = np.random.random() * sum(weights)

        idx = 0
        prob = weights[0]
        logger.debug(self.name + " " + str(weights))
        while prob <= rand and idx < len(weights):
            idx += 1
            try:
                prob += weights[idx]
            except IndexError as e:
                msg = ("Attempted to get index " + str(idx) + " from " +
                       str(repr(weights)) + ": ")
                e.args = tuple([msg + e.args[0]] + list(e.args[1:]))
                raise


        logger_str = "{name} ({cls}) selecting {mtype} (index {idx})"
        logger.info(logger_str.format(
            name=self.name,
            cls=self.__class__.__name__,
            idx=idx,
            mtype=self.movers[idx].name
        ))

        mover = self.movers[idx]

        details = MoveDetails()
        details.inputs = []
        details.choice = idx
        details.chosen_mover = mover
        details.probability = weights[idx] / sum(weights)
        details.weights = weights

        path = paths.RandomChoicePathMoveChange(
            mover.move(globalstate),
            mover=self,
            details=details
        )

        return path

class RandomChoiceMover(SelectionMover):
    """
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
    """

    def __init__(self, movers, weights=None):
        super(RandomChoiceMover, self).__init__(movers)

        if weights is None:
            weights = [1.0] * len(movers)

        self.movers = movers
        self.weights = weights

        initialization_logging(init_log, self,
                               entries=['weights'])

    def _selector(self, globalstate):
        return self.weights

class RandomAllowedChoiceMover(RandomChoiceMover):
    """
    Chooses a random mover from its movers which have existing samples.

    This is different from random choice moves in that this mover only picks
    from sub movers that actually can succeed because they have samples in all
    required input_ensembles

    Attributes
    ----------
    movers : list of PathMover
        the PathMovers to choose from
    weights : list of floats
        the relative weight of each PathMover (does not need to be normalized)
    """

    def _selector(self, globalstate):
        if self.weights is None:
            weights = [1.0] * len(self.movers)
        else:
            weights = list(self.weights) # make a copy

        # this is implemented by setting all weights locally to zero that
        # correspond to movers that will potentially fail since the required
        # input ensembles are not present in the globalstate

        present_ensembles = globalstate.ensembles

        for idx, mover in enumerate(self.movers):
            for ens in mover.input_ensembles:
                if ens not in present_ensembles:
                    # ens might be required but is not present
                    weights[idx] = 0.0

        return weights

class FirstAllowedMover(SelectionMover):
    """
    Chooses a first mover that has samples in all required ensembles.

    A mover can only safely be run, if all inputs can be satisfied. This will pick
    the first mover from the list where all ensembles from input_ensembles are
    found.

    Attributes
    ----------
    movers : list of PathMover
        the PathMovers to choose from
    """

    def _selector(self, globalstate):
        weights = [1.0] * len(self.movers)

        present_ensembles = globalstate.ensembles

        found = False

        for idx, mover in enumerate(self.movers):
            if not found:
                for ens in mover.input_ensembles:
                    if ens not in present_ensembles:
                        # ens might be required but is not present
                        weights[idx] = 0.0

                if weights[idx] > 0.0:
                    found = True
            else:
                weights[idx] = 0.0

        return weights

class LastAllowedMover(SelectionMover):
    """
    Chooses the last mover that has samples in all required ensembles.

    A mover can only safely be run, if all inputs can be satisfied. This will pick
    the last mover from the list where all ensembles from input_ensembles are
    found.

    Attributes
    ----------
    movers : list of PathMover
        the PathMovers to choose from
    """

    def _selector(self, globalstate):
        weights = [1.0] * len(self.movers)

        present_ensembles = globalstate.ensembles

        found = False

        for idx, mover in reversed(list(enumerate(self.movers))):
            if not found:
                for ens in mover.input_ensembles:
                    if ens not in present_ensembles:
                        # ens might be required but is not present
                        weights[idx] = 0.0

                if weights[idx] > 0.0:
                    found = True
            else:
                weights[idx] = 0.0

        return weights




class ConditionalMover(PathMover):
    """
    An if-then-else structure for PathMovers.

    Returns a SequentialPathMoveChange of the if_move movepath and the then_move
    movepath (if if_move is accepted) or the else_move movepath (if if_move
    is rejected).
    """
    def __init__(self, if_mover, then_mover, else_mover):
        """
        Parameters
        ----------
        if_mover : openpathsampling.PathMover
        then_mover : openpathsampling.PathMover
        else_mover : openpathsampling.PathMover
        """
        super(ConditionalMover, self).__init__()
        self.if_mover = if_mover
        self.then_mover = then_mover
        self.else_mover = else_mover
        initialization_logging(init_log, self,
                               ['if_mover', 'then_mover', 'else_mover'])

    @property
    def submovers(self):
        return [self.if_mover, self.then_mover, self.else_mover]

    def _get_in_ensembles(self):
        return [ sub.input_ensembles for sub in self.submovers ]

    def _get_out_ensembles(self):
        return [ sub.output_ensembles for sub in self.submovers ]

    def move(self, globalstate):
        subglobal = globalstate

        ifclause = self.if_mover.move(subglobal)
        samples = ifclause.results
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
    """
    Performs each of the moves in its movers list. Returns all samples
    generated, in the order of the mover list.

    For example, this would be used to create a move that does a sequence of
    replica exchanges in a given order, regardless of whether the moves
    succeed or fail.
    """
    def __init__(self, movers):
        """
        Parameters
        ----------
        movers : list of openpathsampling.PathMover
            the list of pathmovers to be run in sequence
        """
        super(SequentialMover, self).__init__()
        self.movers = movers
        initialization_logging(init_log, self, ['movers'])

    @property
    def submovers(self):
        return self.movers

    @property
    def is_ensemble_change_mover(self):
        if self._is_ensemble_change_mover is not None:
            return self._is_ensemble_change_mover
        sub_change = False
        for mover in self.movers:
            if mover.is_ensemble_change_mover:
                sub_change = True
                break
        return sub_change


    def _get_in_ensembles(self):
        return [ sub.input_ensembles for sub in self.submovers ]

    def _get_out_ensembles(self):
        return [ sub.output_ensembles for sub in self.submovers ]

    def move(self, globalstate):
        logger.debug("Starting sequential move")

        subglobal = globalstate
        pathmovechanges = []

        for mover in self.movers:
            logger.debug("Starting sequential move step "+str(mover))

            # Run the sub mover
            movepath = mover.move(subglobal)
            samples = movepath.results
            subglobal = subglobal.apply_samples(samples)
            pathmovechanges.append(movepath)

        return paths.SequentialPathMoveChange(pathmovechanges, mover=self)


class PartialAcceptanceSequentialMover(SequentialMover):
    """
    Performs each move in its movers list until complete or until one is not
    accepted. If any move is not accepted, further moves are not attempted,
    but the previous accepted samples remain accepted.

    For example, this would be used to create a bootstrap promotion move,
    which starts with a shooting move, followed by an EnsembleHop/Replica
    promotion ConditionalSequentialMover. Even if the EnsembleHop fails, the
    accepted shooting move should be accepted.
    """
    def move(self, globalstate):
        logger.debug("==== BEGINNING " + self.name + " ====")
        subglobal = paths.SampleSet(self.legal_sample_set(globalstate))
        pathmovechanges = []
        for mover in self.movers:
            logger.info(str(self.name)
                        + " starting mover index " + str(self.movers.index(mover) )
                        + " (" + mover.name + ")"
                       )
            # Run the sub mover
            movepath = mover.move(subglobal)
            samples = movepath.results
            subglobal = subglobal.apply_samples(samples)
            pathmovechanges.append(movepath)
            if not movepath.accepted:
                break

        logger.debug("==== FINISHING " + self.name + " ====")
        return paths.PartialAcceptanceSequentialPathMoveChange(pathmovechanges, mover=self)



class ConditionalSequentialMover(SequentialMover):
    """
    Performs each move in its movers list until complete or until one is not
    accepted. If any move in not accepted, all previous samples are updated
    to have set their acceptance to False.

    For example, this would be used to create a minus move, which consists
    of first a replica exchange and then a shooting (extension) move. If the
    replica exchange fails, the move is aborted before doing the dynamics.

    ConditionalSequentialMover only works if there is a *single* active
    sample per replica.
    """
    def move(self, globalstate):
        logger.debug("Starting conditional sequential move")

        subglobal = globalstate
        pathmovechanges = []

        for mover in self.movers:
            logger.debug("Starting sequential move step "+str(mover))

            # Run the sub mover
            movepath = mover.move(subglobal)
            samples = movepath.results
            subglobal = subglobal.apply_samples(samples)
            pathmovechanges.append(movepath)

            if not movepath.accepted:
                break

        return paths.ConditionalSequentialPathMoveChange(pathmovechanges, mover=self)


# TODO: Restrict to last should not be used, but rather a filter by ensemble.
# reason is that the order or samples is partially arbitrary and so the result
# of this mover depends on the implementation of the preceeding mover!!!
# Hence, it might cause hard to find errors!
class RestrictToLastSampleMover(PathMover):
    def __init__(self, mover):
        super(RestrictToLastSampleMover, self).__init__()
        self.mover = mover

    @property
    def submovers(self):
        return [self.mover]

    def _get_in_ensembles(self):
        return [ sub.input_ensembles for sub in self.submovers ]

    def move(self, globalstate):
        movepath = self.mover.move(globalstate)
        return paths.KeepLastSamplePathMoveChange(movepath, mover=self)


class ReplicaIDChangeMover(PathMover):
    """
    Changes the replica ID for a path.
    """
    def __init__(self, replica_pair):
        super(ReplicaIDChangeMover, self).__init__()
        self.replica_pair = replica_pair
        initialization_logging(logger=init_log, obj=self,
                               entries=['replica_pairs'])

    def move(self, globalstate):
        rep_from = self.replica_pair[0]
        rep_to = self.replica_pair[1]
        rep_sample = self.select_sample(globalstate,
                                        ensembles=None,
                                        replicas=rep_from)

        logger.info("Creating new sample from replica ID " + str(rep_from)
                    + " and putting it in replica ID " + str(rep_to))

        # note: currently this clones into a new replica ID. We might later
        # want to kill the old replica ID (and possibly rename this mover).

        sample_details = SampleDetails()

        new_sample = paths.Sample(
            replica=rep_to,
            ensemble=rep_sample.ensemble,
            trajectory=rep_sample.trajectory,
            parent=rep_sample,
            mover=self
        )

        # Can be used to remove the old sample. Not used yet!
        kill_sample = paths.Sample(
            replica=rep_from,
            trajectory=None,
            ensemble=rep_sample.ensemble,
            parent=None,
            mover=self
        )

        details = MoveDetails()
        details.inputs = [rep_sample]
        details.trials = [rep_sample]
        details.mover = self
        setattr(details, 'rep_from', rep_from)
        setattr(details, 'rep_to', rep_to)

        return paths.AcceptedSamplePathMoveChange(
            samples=[new_sample],
            mover=self,
            details=details
        )


class SubPathMover(PathMover):
    """Mover that delegates to a single submover
    """
    def __init__(self, mover):
        """
        Parameters
        ----------
        mover : PathMover
            the submover to be delegated to
        ensembles : nested list of Ensemble or None
            the ensemble specification
        """
        super(SubPathMover, self).__init__()
        self.mover = mover

    @property
    def submovers(self):
        return [self.mover]

    @property
    def is_ensemble_change_mover(self):
        return self.mover.is_ensemble_change_mover

    def _get_in_ensembles(self):
        return self.mover.input_ensembles

    def _get_out_ensembles(self):
        return self.mover.output_ensembles

    def move(self, globalstate):
        subchange = self.mover.move(globalstate)
        change = paths.SubPathMoveChange(
            subchange=subchange,
            mover=self
        )
        return change

#    @classmethod
#    def from_dict(cls, dct):
#        # This will always fix the mover to be the one stored for all SubPathMovers
#        obj = PathMover.from_dict(dct)
#        obj.mover = dct['mover']
#
#        return obj

class EnsembleFilterMover(SubPathMover):
    """Mover that return only samples from specified ensembles
    """
    def __init__(self, mover, ensembles):
        """
        Parameters
        ----------
        mover : PathMover
            the submover to be delegated to
        ensembles : nested list of Ensemble or None
            the ensemble specification
        """
        super(SubPathMover, self).__init__()
        self.ensembles = ensembles
        self.mover = mover


        if not set(self.mover.output_ensembles) & set(self.ensembles):
            # little sanity check, if the underlying move will be removed by the
            # filter throw a warning
            raise ValueError('Your filter removes the underlying move completely. ' +
                             'Please check your ensembles and submovers!')

    def move(self, globalstate):
        # TODO: This will only pass filtered samples. We might split this into an
        # separate input and output filter if only one side is needed

        filtered_globalstate = paths.SampleSet([
            samp for samp in globalstate if samp.ensemble in self.ensembles
        ])
        subchange = self.mover.move(filtered_globalstate)
        change = paths.FilterByEnsemblePathMoveChange(
            subchange=subchange,
            mover=self
        )
        return change

    def _get_in_ensembles(self):
        # only filter the output, not the input
        # return self.mover.input_ensembles
        return self.ensembles

    def _get_out_ensembles(self):
        return self.ensembles


class OneWayShootingMover(RandomChoiceMover):
    """
    OneWayShootingMover is a special case of a RandomChoiceMover which
    combines gives a 50/50 chance of selecting either a ForwardShootMover or
    a BackwardShootMover. Both submovers use the same shooting point
    selector, and both apply to the same ensembles and replicas.

    Attributes
    ----------
    selector : ShootingPointSelector
        The shooting point selection scheme
    ensemble : paths.Ensemble
        Ensemble for this shooting mover
    """
    def __init__(self, ensemble, selector):
        movers = [
            ForwardShootMover(
                ensemble=ensemble,
                selector=selector
            ),
            BackwardShootMover(
                ensemble=ensemble,
                selector=selector
            )
        ]
        super(OneWayShootingMover, self).__init__(
            movers=movers
        )

    @classmethod
    def from_dict(cls, dct):
        mover = cls.__new__(cls)

        # override with stored movers and use the init of the super class
        # this assumes that the super class has movers as its signature
        super(cls, mover).__init__(
            movers=dct['movers']
        )

        return mover

    @property
    def ensemble(self):
        return self.movers[0].ensemble

    @property
    def selector(self):
        return self.movers[0].selector

class OneWayExtendMover(RandomChoiceMover):
    """
    OneWayShootingMover is a special case of a RandomChoiceMover which
     gives a 50/50 chance of selecting either a ForwardExtendMover or
    a BackwardExtendMover. Both submovers use the same same ensembles
    and replicas.

    Attributes
    ----------
    ensembles : openpathsampling.Ensemble
        valid ensemble
    """
    def __init__(self, ensemble, target_ensemble):
        movers = [
            ForwardExtendMover(
                ensemble=ensemble,
                target_ensemble=target_ensemble
            ),
            BackwardExtendMover(
                ensemble=ensemble,
                target_ensemble=target_ensemble
            )
        ]
        super(OneWayExtendMover, self).__init__(
            movers=movers
        )

    @classmethod
    def from_dict(cls, dct):
        mover = cls.__new__(cls)

        # override with stored movers and use the init of the super class
        # this assumes that the super class has movers as its signature
        super(cls, mover).__init__(
            movers=dct['movers']
        )

        return mover

class MinusMover(SubPathMover):
    """
    Instance of a MinusMover.

    The minus move combines a replica exchange with path extension to swap
    paths between the innermost regular TIS interface ensemble and the minus
    interface ensemble. This is particularly useful for improving sampling
    of path space.
    """
    _is_canonical = True

    def __init__(self, minus_ensemble, innermost_ensembles):
        try:
            innermost_ensembles = list(innermost_ensembles)
        except TypeError:
            innermost_ensembles = [innermost_ensembles]

        segment = minus_ensemble._segment_ensemble
        sub_trajectory_selector = RandomChoiceMover([
            FirstSubtrajectorySelectMover(
                ensemble=minus_ensemble,
                sub_ensemble=segment,
                n_l=minus_ensemble.n_l
            ),
            FinalSubtrajectorySelectMover(
                ensemble=minus_ensemble,
                sub_ensemble=segment,
                n_l=minus_ensemble.n_l
                ),
        ])
        sub_trajectory_selector.name = "MinusSubtrajectoryChooser"

        repexs = [ReplicaExchangeMover(
            ensemble1=segment,
            ensemble2=inner
        ) for inner in innermost_ensembles]

        repex_chooser = RandomChoiceMover(repexs)
        repex_chooser.name = "InterfaceSetChooser"

        extension_mover = RandomChoiceMover([
            ForwardExtendMover(
                ensemble=segment,
                target_ensemble=minus_ensemble
            ),
            BackwardExtendMover(
                ensemble=segment,
                target_ensemble=minus_ensemble
            )
        ])

        extension_mover.name = "MinusExtensionDirectionChooser"
        self.engine = extension_mover.movers[0].engine
        if self.engine is not extension_mover.movers[1].engine:
            raise RuntimeWarning("Forward and backward engines differ?!?!")

        mover = \
            EnsembleFilterMover(
                ConditionalSequentialMover([
                    sub_trajectory_selector,
                    repex_chooser,
                    extension_mover
                ]),
            ensembles=[minus_ensemble] + innermost_ensembles
        )

        self.minus_ensemble = minus_ensemble
        self.innermost_ensembles = innermost_ensembles
        initialization_logging(init_log, self, ['minus_ensemble',
                                                'innermost_ensembles'])

        super(MinusMover, self).__init__(mover)

class SingleReplicaMinusMover(MinusMover):
    """
    Minus mover for single replica TIS.

    In SRTIS, the minus mover doesn't actually keep an active sample in the
    minus interface. Instead, it just puts the newly generated segment into
    the innermost ensemble.
    """
    def __init__(self, minus_ensemble, innermost_ensembles, bias=None):
        try:
            innermost_ensembles = list(innermost_ensembles)
        except TypeError:
            innermost_ensembles = [innermost_ensembles]

        # TODO: Until we have automated detailed balance calculations, I
        # think this will only be valid in the case of only one innermost
        # ensemble.  But I think you only want to use it in the case of only
        # one innermost ensemble anyway. The following warns us:
        if len(innermost_ensembles) > 1:
            logger.warning("Probably shouldn't use SingleReplicaMinusMover with MISTIS")

        segment = minus_ensemble._segment_ensemble

        hop_innermost_to_segment = RandomAllowedChoiceMover([
            EnsembleHopMover(innermost, segment, bias=bias)
            for innermost in innermost_ensembles
        ])

        # TODO: again, works for single interface set, but there has to be a
        # smarter way to do this in the MISTIS case
        hop_segment_to_innermost = RandomChoiceMover([
            EnsembleHopMover(segment, innermost, bias=bias)
            for innermost in innermost_ensembles
        ])

        forward_minus = ConditionalSequentialMover([
            hop_innermost_to_segment,
            ForwardExtendMover(segment, minus_ensemble),
            FinalSubtrajectorySelectMover(minus_ensemble, segment),
            hop_segment_to_innermost
        ])

        backward_minus = ConditionalSequentialMover([
            hop_innermost_to_segment,
            BackwardExtendMover(segment, minus_ensemble),
            FirstSubtrajectorySelectMover(minus_ensemble, segment),
            hop_segment_to_innermost
        ])

        mover = EnsembleFilterMover(RandomChoiceMover([backward_minus, 
                                                       forward_minus]),
                                    ensembles=innermost_ensembles)

        # we skip MinusMover's init and go to the grandparent
        super(MinusMover, self).__init__(mover)



class PathSimulatorMover(SubPathMover):
    """
    This just wraps a mover and references the used pathsimulator
    """
    def __init__(self, mover, pathsimulator):
        super(PathSimulatorMover, self).__init__(mover)
        self.pathsimulator = pathsimulator

    def move(self, globalstate, step=-1):
        details = MoveDetails(
            step=step
        )

        return paths.PathSimulatorPathMoveChange(
            self.mover.move(globalstate),
            mover=self,
            details=details
        )


class MultipleSetMinusMover(RandomChoiceMover):
    pass

def NeighborEnsembleReplicaExchange(ensemble_list):
    movers = [
        ReplicaExchangeMover(
            ensemble1=ensemble_list[i],
            ensemble2=ensemble_list[i+1]
        )
        for i in range(len(ensemble_list)-1)
    ]
    return movers

def PathReversalSet(ensembles):
    return map(PathReversalMover, ensembles)


class PathMoverFactory(object):
    @staticmethod
    def OneWayShootingSet(selector_set, interface_set):
        if type(selector_set) is not list:
            selector_set = [selector_set]*len(interface_set)

        mover_set = []
        for (selector, iface) in zip(selector_set, interface_set):
            mover = OneWayShootingMover(
                selector=selector,
                ensemble=iface
            )
            mover.name = "OneWayShootingMover " + str(iface.name)
            mover_set.append(mover)

        return mover_set

    @staticmethod
    def TwoWayShootingSet():
        pass

    @staticmethod
    def NearestNeighborRepExSet():
        pass


class Details(StorableObject):
    """Details of an object. Can contain any data
    """

    def __init__(self, **kwargs):
        super(Details, self).__init__()
        for key, value in kwargs.iteritems():
            setattr(self, key, value)

    def __str__(self):
        # primarily for debugging/interactive use
        mystr = ""
        for key in self.__dict__.keys():
            if not isinstance(self.__dict__[key], paths.Ensemble):
                mystr += str(key) + " = " + str(self.__dict__[key]) + '\n'
        return mystr


class MoveDetails(Details):
    """Details of the move as applied to a given replica

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
    """

    def __init__(self, **kwargs):
        self.inputs=None
        self.trials=None
        self.results=None
        super(MoveDetails, self).__init__(**kwargs)


class SampleDetails(Details):
    """Details of a sample

    Attributes
    ----------
    selection_probability : float
        the chance that a sample will be accepted due to asymmetrical proposal
    """

    def __init__(self, **kwargs):
        self.selection_probability=1.0
        super(SampleDetails, self).__init__(**kwargs)
