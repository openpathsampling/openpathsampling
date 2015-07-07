"""
Created on 19.07.2014

@author: Jan-Hendrik Prinz
@author: David W. H. Swenson
"""

import numpy as np
import random

import openpathsampling as paths
from openpathsampling.todict import OPSNamed, OPSObject

import logging
from ops_logging import initialization_logging

from treelogic import TreeMixin

logger = logging.getLogger(__name__)
init_log = logging.getLogger('openpathsampling.initialization')

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

    if list_of_lists == True:
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


class PathMover(TreeMixin, OPSNamed):
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


    Attributes
    ----------
    name : string
        a human-readable name of the PathMover instance
    ensembles : nested list of Ensemble or None
        a mover-specific representation of the ensembles the mover acts on.
        Usually this is either None (meaning all ensembles) or a specific
        set of ensembles

    """

    def __init__(self,  ensembles=None):
        OPSNamed.__init__(self)

        # we keep ensembles totally arbitrary. Each mover know what to do
        self.ensembles = ensembles

        self._in_ensembles = None
        self._out_ensembles = None
        self._len = None

        initialization_logging(logger=init_log, obj=self,
                               entries=['ensembles'])

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
        return self.ensembles

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

class SwappingMover(MoverType):
    """
    A mover that swaps samples from ensembles in some way. Relevant for mixing
    """


###############################################################################
# GENERATORS
###############################################################################

class SampleGeneratingMover(PathMover):
    engine = None

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
        # GeneratingMovers do not have submovers!
        return []

    def _ensemble_selector(self, globalstate):
        """Function to determine which ensembles to pick samples from

        Returns
        -------
        list of Ensemble
            the list of ensembles. Samples can then be selected using
            PathMover.select_sample
        """

        # Default is that the list of ensembles is in self.ensembles
        return self.ensembles

    def __init__(self, ensembles=None):
        super(SampleGeneratingMover, self).__init__(ensembles)

    def move(self, globalstate):
        # 1. pick a set of ensembles (in case we allow to pick several ones)
        ensembles = self._ensemble_selector(globalstate)

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
        then used to generate new samples. If the GeneratingMover is used as a move
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

class EngineGeneratingMover(SampleGeneratingMover):
    """Baseclass for GeneratingMovers that use an engine
    """

    engine = None

class ShootGeneratingMover(EngineGeneratingMover):
    """Main class for GeneratingMovers using ShootingMoves

    Attributes
    ----------
    selector
    """
    def __init__(self, selector, ensembles=None):
        """
        Parameters
        ----------
        selector : ShootingPointSelector
            the shootingpoint selector to determine the shooting point in the
            move
        ensembles : Ensemble of None
            the specific ensemble to be shot from or None if all are allowed
        """
        super(ShootGeneratingMover, self).__init__(ensembles)
        self.selector = selector

    def _ensemble_selector(self, globalstate):
        # return a single ensemble
        return [ self.ensembles ]

    def __call__(self, trial):
        initial_trajectory = trial.trajectory

        dynamics_ensemble = trial.ensemble
        replica = trial.replica

        initial_point = self.selector.pick(initial_trajectory)
        trial_point = self._shoot(initial_point, dynamics_ensemble)

        bias = initial_point.sum_bias / trial_point.sum_bias

        trial_details = paths.SampleDetails(
            initial_point=initial_point,
            trial_point=trial_point,
        )

        trial = paths.Sample(
            replica=replica,
            trajectory=trial_point.trajectory,
            ensemble=dynamics_ensemble,
            parent=trial,
            details=trial_details,
            mover=self,
            bias=bias
        )

        trials = [trial]

        return trials

    def _shoot(self, shooting_point, ensemble):
        """Implementation of the shooting

        Parameters
        ----------
        shooting_point : ShootingPoint
            the initial shooting point instance containing a reference to the
            initial trajectory
        ensemble : Ensemble
            the ensemble for which the trajectory is to be created. This defines
            the stopping criterion

        Returns
        -------
        ShootingPoint
            the final shooting point referencing the final trajectory
        """
        return shooting_point


class ForwardShootGeneratingMover(ShootGeneratingMover):
    """A forward shooting sample generator
    """
    def _shoot(self, shooting_point, ensemble):
        shoot_str = "Shooting {sh_dir} from frame {fnum} in [0:{maxt}]"
        logger.info(shoot_str.format(fnum=shooting_point.index,
                                     maxt=len(shooting_point.trajectory)-1,
                                     sh_dir="forward",
                                    ))

        # Run until one of the stoppers is triggered
        partial_trajectory = self.engine.generate(
            shooting_point.snapshot.copy(),
            running = [
                paths.PrefixTrajectoryEnsemble(
                    ensemble,
                    shooting_point.trajectory[0:shooting_point.index]
                ).can_append,
                self.engine.max_length_stopper.can_append
            ]
        )

        trial_trajectory = shooting_point.trajectory[0:shooting_point.index] + partial_trajectory

        trial_point = paths.ShootingPoint(
            shooting_point.selector,
            trial_trajectory,
            shooting_point.index
        )

        return trial_point


class BackwardShootGeneratingMover(ShootGeneratingMover):
    """A Backward shooting generator
    """
    def _shoot(self, shooting_point, ensemble):
        shoot_str = "Shooting {sh_dir} from frame {fnum} in [0:{maxt}]"
        logger.info(shoot_str.format(fnum=shooting_point.index,
                                     maxt=len(shooting_point.trajectory)-1,
                                     sh_dir="backward",
                                    ))

        # Run until one of the stoppers is triggered
        partial_trajectory = self.engine.generate(
            shooting_point.snapshot.reversed_copy(),
            running = [
                paths.SuffixTrajectoryEnsemble(
                    ensemble,
                    shooting_point.trajectory[shooting_point.index + 1:]
                ).can_prepend,
                self.engine.max_length_stopper.can_append
            ]
        )

        trial_trajectory = partial_trajectory.reversed + shooting_point.trajectory[shooting_point.index + 1:]

        trial_point = paths.ShootingPoint(
            shooting_point.selector,
            trial_trajectory,
            len(partial_trajectory) - 1
        )

        return trial_point

# TODO: This doubling might be superfluous


class ShootMover(ShootGeneratingMover):
    """
    A pathmover that implements a general shooting algorithm
    """

class ForwardShootMover(ForwardShootGeneratingMover):
    """
    A pathmover that implements the forward shooting algorithm
    """


class BackwardShootMover(BackwardShootGeneratingMover):
    """
    A pathmover that implements the backward shooting algorithm
    """


###############################################################################
# EXTENDING GENERATORS
###############################################################################

class ExtendingGeneratingMover(EngineGeneratingMover):
    """
    Sample GeneratingMover that creates Samples using extensions

    Extending will create samples in a super ensemble from samples
    in a smaller ensemble by forward or backward extending the original
    sample until it is in the target ensemble. This requires the the target
    ensemble is reachable from the initial ensemble
    """
    def __init__(self, extend_ensemble, ensembles=None):
        """
        Parameters
        ----------
        extend_ensemble : openpathsampling.Ensemble
            the target ensemble
        ensembles : openpathsampling.Ensemble
            the initial ensemble to be started from
        """
        super(ExtendingGeneratingMover, self).__init__(
            ensembles
        )
        self.extend_ensemble = extend_ensemble

    def _ensemble_selector(self, globalstate):
        return [self.ensembles]

    def _get_in_ensembles(self):
        return [self.ensembles]

    def _get_out_ensembles(self):
        return [self.extend_ensemble]

    def __call__(self, trial):
        initial_trajectory = trial.trajectory

        replica = trial.replica
        trial_trajectory = self._extend(initial_trajectory, self.extend_ensemble)

        trial_details = paths.SampleDetails(
        )

        # the actual bias would be 0.0 since we will never be able to do the
        # reverse move. Since this is the opposite of subtraj we set both
        # proposal bias for these to 100% which means no bias

        trial = paths.Sample(
            replica=replica,
            trajectory=trial_trajectory,
            ensemble=self.extend_ensemble,
            parent=trial,
            details=trial_details,
            mover=self,
            bias=1.0
        )

        trials = [trial]

        return trials

    def _extend(self, initial_trajectory, ensemble):
        return initial_trajectory


class ForwardExtendGeneratingMover(ExtendingGeneratingMover):
    """
    A Sample GeneratingMover implementing Forward Extension
    """
    def _extend(self, initial_trajectory, ensemble):
        shoot_str = "Extending {sh_dir} from frame {fnum} in [0:{maxt}]"
        logger.info(shoot_str.format(fnum=len(initial_trajectory)-1,
                                     maxt=len(initial_trajectory)-1,
                                     sh_dir="forward",
                                    ))

        # Run until one of the stoppers is triggered
        partial_trajectory = self.engine.generate(
            initial_trajectory[-1],
            running = [
                paths.PrefixTrajectoryEnsemble(
                    ensemble,
                    initial_trajectory[:-1]
                ).can_append,
                self.engine.max_length_stopper.can_append
            ]
        )

        trial_trajectory = initial_trajectory + partial_trajectory[1:]

        return trial_trajectory


class BackwardExtendGeneratingMover(ExtendingGeneratingMover):
    """
    A Sample GeneratingMover implementing Backward Extension
    """
    def _extend(self, initial_trajectory, ensemble):
        shoot_str = "Extending {sh_dir} from frame {fnum} in [0:{maxt}]"
        logger.info(shoot_str.format(fnum=0,
                                     maxt=len(initial_trajectory)-1,
                                     sh_dir="backward",
                                    ))

        # Run until one of the stoppers is triggered
        partial_trajectory = self.engine.generate(
            initial_trajectory[0].reversed,
            running = [
                paths.SuffixTrajectoryEnsemble(
                    ensemble,
                    initial_trajectory[1:]
                ).can_prepend,
                self.engine.max_length_stopper.can_append
            ]
        )

        trial_trajectory = partial_trajectory.reversed + initial_trajectory[1:]
        return trial_trajectory


class ForwardExtendMover(ForwardExtendGeneratingMover):
    """Creates a new sample by extending forward to a new ensemble
    """


class BackwardExtendMover(BackwardExtendGeneratingMover):
    """Creates a new sample by extending backward to a new ensemble
    """


###############################################################################
# REPLICA EXCHANGE GENERATORS
###############################################################################

class ReplicaExchangeGeneratingMover(SampleGeneratingMover):
    """
    A Sample GeneratingMover implementing a standard Replica Exchange
    """
    _is_ensemble_change_mover = True

    def __init__(self, bias=None, ensembles=None):
        """
        Parameters
        ----------
        bias : list of float
            bias is not used yet
        ensembles : list of openpathsampling.Ensemble
        """
        # either replicas or ensembles must be a list of pairs; more
        # complicated filtering can be done with a wrapper class
        super(ReplicaExchangeGeneratingMover, self).__init__(ensembles)
        # TODO: add support for bias; cf EnsembleHopMover
        self.bias = bias
        initialization_logging(logger=init_log, obj=self,
                               entries=['bias'])

    def _ensemble_selector(self, globalstate):
        list_of_ensemble_pairs = make_list_of_pairs(self.ensembles)
        selected = random.choice(list_of_ensemble_pairs)
        return selected

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

class StateSwapGeneratingMover(SampleGeneratingMover):
    def __init__(self, bias=None, ensembles=None):
        # either replicas or ensembles must be a list of pairs; more
        # complicated filtering can be done with a wrapper class
        super(StateSwapGeneratingMover, self).__init__(ensembles)
        self.bias = bias
        initialization_logging(logger=init_log, obj=self,
                               entries=['bias'])

    def _ensemble_selector(self, globalstate):
        list_of_ensemble_pairs = make_list_of_pairs(self.ensembles)
        selected = random.choice(list_of_ensemble_pairs)
        return selected

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

class ReplicaExchangeMover(ReplicaExchangeGeneratingMover):
    pass

class StateSwapMover(StateSwapGeneratingMover):
    pass

###############################################################################
# SUBTRAJECTORY GENERATORS
###############################################################################


class RandomSubtrajectorySelectGeneratingMover(SampleGeneratingMover):
    """
    Samples a random subtrajectory satisfying the given subensemble.

    If there are no subtrajectories which satisfy the subensemble, this
    returns the zero-length trajectory.

    Parameters
    ----------
    subensemble : Ensemble
        the subensemble to be searched for
    n_l : int
        the number of
    ensembles : list of Ensembles or None
        the set of allows samples to chose from

    Attributes
    ----------


    """
    _is_ensemble_change_mover = True
    def __init__(self, subensemble, n_l=None, ensembles=None):
        super(RandomSubtrajectorySelectGeneratingMover, self).__init__(
            ensembles
        )
        self.n_l = n_l
        self.subensemble = subensemble

    def _ensemble_selector(self, globalstate):
        return [ self.ensembles ]

    def _get_in_ensembles(self):
        return [ self.ensembles ]

    def _get_out_ensembles(self):
        return [ self.subensemble ]

    def _choose(self, trajectory_list):
        return random.choice(trajectory_list)

    def __call__(self, trial):
        initial_trajectory = trial.trajectory
        replica = trial.replica
        logger.debug("Working with replica " + str(replica) + " (" + str(initial_trajectory) + ")")

        subtrajs = self.subensemble.split(initial_trajectory)
        logger.debug("Found "+str(len(subtrajs))+" subtrajectories.")

        if (self.n_l is None and len(subtrajs) > 0) or \
            (self.n_l is not None and len(subtrajs) == self.n_l):
            subtraj = self._choose(subtrajs)

            bias = 1.0

            trial = paths.Sample(
                replica=replica,
                trajectory=subtraj,
                ensemble=self.subensemble,
                parent=trial,
                mover=self,
                bias=bias
            )

            trials = [trial]
        else:
            trials = []

        return trials


class RandomSubtrajectorySelectMover(RandomSubtrajectorySelectGeneratingMover):
    """
    Samples a random subtrajectory satisfying the given subensemble.

    If there are no subtrajectories which satisfy the subensemble, this
    returns the zero-length trajectory.
    """


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

class PathReversalGeneratingMover(SampleGeneratingMover):

    def _ensemble_selector(self, globalstate):
        return [ self.ensembles ]

    def _get_in_ensembles(self):
        return [ self.ensembles ]

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


class PathReversalMover(PathReversalGeneratingMover):
    pass


class EnsembleHopGeneratingMover(SampleGeneratingMover):
    _is_ensemble_change_mover = True
    def __init__(self, bias=None, ensembles=None):
        """
        Parameters
        ----------
        bias : float, dict or None (default)
            gives the bias of accepting (not proposing) a hop. A float will
            be the acceptance for all possible attempts. If a dict is given,
            then it contains a list of ensembles and a matrix. None means
            no bias
        ensembles : list of ensemble pairs
            the list of possible hop attempts

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
        # TODO: maybe allow a version of this with a single ensemble and ANY
        # ensemble can hop to that? messy to code; maybe same idea under
        # another name
        ensembles = make_list_of_pairs(ensembles)
        super(EnsembleHopGeneratingMover, self).__init__(ensembles=ensembles)
        # ensembles -- another version might take a value for each ensemble,
        # and use the ratio; this latter is better for CITIS
        self.bias = bias
        initialization_logging(
            logger=init_log,
            obj=self,
            entries=['bias']
        )

    def _ensemble_selector(self, globalstate):
        # Picks a random initial ensemble from all possible ones
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
        return [ random.choice(legal_ensembles) ]

    @property
    def submovers(self):
        return []

    def _get_in_ensembles(self):
        return [pair[0] for pair in self.ensembles]

    def _get_out_ensembles(self):
        return [pair[1] for pair in self.ensembles]

    def __call__(self, rep_sample):
        ens_from = rep_sample.ensemble

        # pick a random hop to an allowed final ensemble
        legal_pairs = [pair for pair in self.ensembles
                       if pair[0] is ens_from]
        logger.debug("Legal pairs: " + str(legal_pairs))
        ens_pair = random.choice(legal_pairs)
        ens_to = ens_pair[1]

        logger.debug("Selected sample: " + repr(rep_sample))
        replica = rep_sample.replica

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


class EnsembleHopMover(EnsembleHopGeneratingMover):
    """
    A Mover describing a trial change of ensembles
    """


class RandomChoiceMover(PathMover):
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
    def __init__(self, movers, ensembles=None,  weights=None, name=None):
        super(RandomChoiceMover, self).__init__(ensembles=ensembles)

        if name is not None:
            self.name = name

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
        details.probability = self.weights[idx] / sum(self.weights)

        path = paths.RandomChoicePathMoveChange(
            mover.move(globalstate),
            mover=self,
            details=details
        )

        return path


class EnsembleDictionaryMover(PathMover):
    """
    Selects random sample and picks which move to do based its ensemble.


    """
    def __init__(self, ensemble_to_mover_dict):
        ensembles = ensemble_to_mover_dict.keys()
        super(EnsembleDictionaryMover, self).__init__(ensembles=ensembles)

        self.movers = ensemble_to_mover_dict.values()
        self.ensemble_to_mover_dict = ensemble_to_mover_dict

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
        samp = random.choice(globalstate.samples)
        submove = self.ensemble_to_mover_dict[samp.ensemble]
        subset = paths.SampleSet([samp])
        logger_str = "{name} (EnsembleDictionaryMover) selecting {m} (ensemble {ens})"
        logger.info(logger_str.format(name=self.name, m=submove.name,
                                      ens=samp.ensemble.name))
        
        details = MoveDetails()
        details.inputs = []
        details.choice = samp
        details.chosen_mover = submove
        details.probability = 1.0/len(globalstate)

        subchange = submove.move(subset)

        change = paths.SubPathMoveChange(
            subchange=subchange,
            mover=self,
            details=details
        )
        return change



class ConditionalMover(PathMover):
    """
    An if-then-else structure for PathMovers.

    Returns a SequentialPathMoveChange of the if_move movepath and the then_move
    movepath (if if_move is accepted) or the else_move movepath (if if_move
    is rejected).
    """
    def __init__(self, if_mover, then_mover, else_mover, ensembles=None):
        super(ConditionalMover, self).__init__(ensembles=ensembles)
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
    def __init__(self, movers, ensembles=None):
        super(SequentialMover, self).__init__(ensembles=ensembles)
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
    def __init__(self, replica_pairs, ensembles=None):
        self.replica_pairs = make_list_of_pairs(replica_pairs)
        super(ReplicaIDChangeMover, self).__init__(ensembles=ensembles)
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
        setattr(details, 'rep_from', mypair[0])
        setattr(details, 'rep_to', mypair[1])

        return paths.AcceptedSamplePathMoveChange(
            samples=[new_sample],
            mover=self,
            details=details
        )

# TODO: Filter moves are not used at all, do we need these?
# TODO: Turn Filter into real mover with own movechange ?

class FilterByReplica(PathMover):

    def __init__(self, mover, replicas):
        super(FilterByReplica, self).__init__()
        if type(replicas) is not list:
            replicas = [replicas]
        self.replicas = replicas
        self.mover = mover
        # TODO: clean this up
        pass

    def move(self, globalstate):
        filtered_gs = paths.SampleSet(
            [s for s in globalstate if s.replica in self.replicas]
        )
        return self.mover.move(filtered_gs)


class FilterBySample(PathMover):
    def __init__(self, mover, selected_samples):
        super(FilterBySample, self).__init__()
        if type(selected_samples) is not list:
            selected_samples = [selected_samples]
        self.selected_samples = selected_samples
        self.mover = mover

    def move(self, globalstate):
        return paths.FilterSamplesPathMoveChange(
            self.mover.move(globalstate),
            mover=self
        )


class WrappedMover(PathMover):
    """Mover that delegates to a single submover
    """
    def __init__(self, mover, ensembles=None):
        """
        Parameters
        ----------
        mover : PathMover
            the submover to be delegated to
        ensembles : nested list of Ensemble or None
            the ensemble specification
        """
        super(WrappedMover, self).__init__(ensembles)
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


class EnsembleFilterMover(WrappedMover):
    """Mover that return only samples from specified ensembles
    """

    def move(self, globalstate):
        subchange = self.mover.move(globalstate)
        change = paths.FilterByEnsemblePathMoveChange(
            subchange=subchange,
            mover=self
        )
        return change

    def _get_in_ensembles(self):
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
    ensembles : list of Ensemble or None
        valid ensembles; None implies all ensembles are allowed (no
        restriction)
    """
    def __init__(self, selector, ensembles=None):
        movers = [
            ForwardShootMover(selector, ensembles),
            BackwardShootMover(selector, ensembles)
        ]
        super(OneWayShootingMover, self).__init__(
            movers=movers, ensembles=ensembles
        )
        self.selector = selector


class MinusMover(WrappedMover):
    """
    Instance of a MinusMover.

    The minus move combines a replica exchange with path extension to swap
    paths between the innermost regular TIS interface ensemble and the minus
    interface ensemble. This is particularly useful for improving sampling
    of path space.
    """
    _is_canonical = True

    def __init__(self, minus_ensemble, innermost_ensembles, ensembles=None):
        segment = minus_ensemble._segment_ensemble
        try:
            innermost_ensembles = list(innermost_ensembles)
        except TypeError:
            innermost_ensembles = [innermost_ensembles]
        innermost_ensemble = paths.join_ensembles(innermost_ensembles)
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

        repexs = [ReplicaExchangeMover(ensembles=[[segment, inner]])
                  for inner in innermost_ensembles]
        repex_chooser = RandomChoiceMover(repexs)
        repex_chooser.name = "InterfaceSetChooser"

        extension_mover = RandomChoiceMover([
            ForwardExtendMover(minus_ensemble, segment),
            BackwardExtendMover(minus_ensemble, segment)
        ])

        extension_mover.name = "MinusExtensionDirectionChooser"
        self.engine = extension_mover.movers[0].engine
        if self.engine is not extension_mover.movers[1].engine:
            raise RuntimeWarning("Forward and backward engines differ?!?!")

        mover = EnsembleFilterMover(
            ConditionalSequentialMover([
                subtrajectory_selector,
                repex_chooser,
                extension_mover
            ]),
            ensembles=[minus_ensemble] + innermost_ensembles
        )

        self.minus_ensemble = minus_ensemble
        self.innermost_ensemble = innermost_ensemble
        self.innermost_ensembles = innermost_ensembles
        initialization_logging(init_log, self, ['minus_ensemble',
                                                'innermost_ensemble'])

        super(MinusMover, self).__init__(mover)


class PathSimulatorMover(WrappedMover):
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
        ReplicaExchangeMover(ensembles=[[ensemble_list[i], ensemble_list[i+1]]])
        for i in range(len(ensemble_list)-1)
    ]
    return movers

def PathReversalSet(l):
    # TODO: Check if replica can be removed here
#    if isinstance(l[0], paths.Ensemble):
    return [PathReversalMover(ensembles=[item]) for item in l]
#    else:
#        return [PathReversalMover(replicas=[item]) for item in l]


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
    """Details of an object. Can contain any data
    """

    def __init__(self, **kwargs):
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
