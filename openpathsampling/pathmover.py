"""
Created on 19.07.2014

@author: Jan-Hendrik Prinz
@author: David W. H. Swenson
"""
import abc
import logging
import numpy as np
import random

import openpathsampling as paths
from openpathsampling.netcdfplus import StorableNamedObject, StorableObject
from openpathsampling.pathmover_inout import InOutSet, InOut
from openpathsampling.rng import default_rng
from .ops_logging import initialization_logging
from .treelogic import TreeMixin

from openpathsampling.deprecations import deprecate, has_deprecations
from openpathsampling.deprecations import SAMPLE_DETAILS, MOVE_DETAILS

from future.utils import with_metaclass

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

    Parameters
    ----------
    l : list
        input list, either flat list of length 2N, a list of pairs or None

    Returns
    -------
    list of pairs
    """
    if l is None:
        return None

    _ = len(l)  # raises TypeError, avoids everything else

    # based on first element, decide whether this should be a list of lists
    # or a flat list
    try:
        _ = len(l[0])
        list_of_lists = True
    except TypeError:
        list_of_lists = False

    if list_of_lists:
        for elem in l:
            assert len(elem) == 2, "List of lists: inner list length != 2"
        outlist = l
    else:
        assert len(l) % 2 == 0, "Flattened list: length not divisible by 2"
        outlist = [
            [a, b] for (a, b) in zip(l[slice(0, None, 2)], l[slice(1, None, 2)])
        ]
    # Note that one thing we don't check is whether the items are of the
    # same type. That might be worth doing someday; for now, we trust that
    # part to work.
    return outlist


class SampleNaNError(Exception):
    def __init__(self, message, trial_sample, details):
        super(SampleNaNError, self).__init__(message)
        self.trial_sample = trial_sample
        self.details = details


class SampleMaxLengthError(Exception):
    def __init__(self, message, trial_sample, details):
        super(SampleMaxLengthError, self).__init__(message)
        self.trial_sample = trial_sample
        self.details = details


class MoveChangeNaNError(Exception):
    pass


class PathMover(with_metaclass(abc.ABCMeta, TreeMixin, StorableNamedObject)):
    """
    A PathMover is the description of a move in replica space.

    Notes
    -----
    A pathmover takes a SampleSet() and returns MoveChange() that is
    used to change the old SampleSet() to the new one.

    SampleSet1 + MoveChange1 => SampleSet2

    A MoveChange is effectively a list of Samples. The change acts upon
    a SampleSet by replacing existing Samples in the same ensemble
    sequentially.

    SampleSet({samp1(ens1), samp2(ens2), samp3(ens3)}) +
        MoveChange([samp4(ens2)])
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

    #__metaclass__ = abc.ABCMeta

    def __init__(self):
        StorableNamedObject.__init__(self)
        self._rng = default_rng()
        self._in_ensembles = None
        self._out_ensembles = None
        self._len = None
        self._inout = None
        self._trust_candidate = False

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

    # +-------------------------------------------------------------------------
    # | analyze effects of sample sets
    # +-------------------------------------------------------------------------

    def move_replica_state(self, replica_states):
        return self.in_out.move(replica_states)

    def sub_replica_state(self, replica_states):
        """
        Return set of replica states that a submover might be called with

        Parameters
        ----------
        replica_states : set of `openpathsampling.pathmover_inout.ReplicaState`

        Returns
        -------
        list of set of `ReplicaState`

        """
        return [replica_states] * len(self.submovers)

    def _generate_in_out(self):
        if len(self.output_ensembles) == 0:
            return {
                InOutSet([])
            }
        elif len(self.input_ensembles) == 1 and len(self.output_ensembles) == 1:
            return InOutSet([
                InOut(
                    [((self.input_ensembles[0], self.output_ensembles[0], 0), 1)]
                )])
        else:
            # Fallback could be all possibilities, but for now we ask the user!
            raise NotImplementedError(
                'Please implement the in-out-matrix for this mover.')

    @property
    def in_out(self):
        """
        List the input -> output relation for ensembles

        A mover will pick one or more replicas from specific ensembles.
        Alter them (or not) and place these (or additional ones) in specific
        ensembles. This relation can be visualized as a mapping of input to
        output ensembles. Like

        ReplicaExchange
        ens1 -> ens2
        ens2 -> ens1

        EnsembleHop (A sample in ens1 will disappear and appear in ens2)
        ens1 -> ens2

        DuplicateMover (create a copy with a new replica number) Not used yet!
        ens1 -> ens1
        None -> ens1

        Returns
        -------
        list of list of tuple : (:obj:`openpathsampling.Ensemble`,
        :obj:`openpathsampling.Ensemble`)
            a list of possible lists of tuples of ensembles.

        Notes
        -----
        The default implementation will
        (1) in case of a single input and output connect the two,
        (2) return nothing if there are no out_ensembles and
        (3) for more then two require implementation
        """
        if self._inout is None:
            self._inout = self._generate_in_out()

        return self._inout

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
        return inp, out

    @property
    def ensemble_signature(self):
        return self._ensemble_signature()

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
        list of :class:`openpathsampling.Ensemble`
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

    @staticmethod
    def legal_sample_set(sample_set, ensembles=None, replicas='all'):
        """
        This returns all the samples from sample_set which are in both
        self.replicas and the parameter ensembles. If ensembles is None, we
        use self.ensembles. If you want all ensembles allowed, pass
        ensembles='all'.

        Parameters
        ----------
        sample_set : `openpathsampling.SampleSet`
            the sampleset from which to pick specific samples matching certain
            criteria
        ensembles : list of `openpathsampling.Ensembles`
            the ensembles to pick from
        replicas : list of int or `all`
            the replicas to pick or `'all'` for all
        """
        mover_replicas = sample_set.replica_list()

        if replicas == 'all':
            selected_replicas = sample_set.replica_list()
        else:
            selected_replicas = replicas

        reps = list(set(mover_replicas) & set(selected_replicas))
        rep_samples = []
        for rep in reps:
            rep_samples.extend(sample_set.all_from_replica(rep))

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
                #     ens_samples.extend(sample_set.all_from_ensemble(ens[0]))
                # except TypeError:
                ens_samples.extend(sample_set.all_from_ensemble(ens))
            legal_samples = list(set(rep_samples) & set(ens_samples))

        return legal_samples

    @staticmethod
    def select_sample(sample_set, ensembles=None, replicas=None):
        """
        Returns one of the legal samples given self.replica and the ensemble
        set in ensembles.

        Parameters
        ----------
        sample_set : `openpathsampling.SampleSet`
            the sampleset from which to pick specific samples matching certain
            criteria
        ensembles : list of `openpathsampling.Ensembles` or `None`
            the ensembles to pick from or `None` for all
        replicas : list of int or None
            the replicas to pick or `None` for all

        """
        if replicas is None:
            replicas = 'all'

        logger.debug(
            "replicas: " + str(replicas) + " ensembles: " + repr(ensembles))
        legal = PathMover.legal_sample_set(sample_set, ensembles, replicas)
        for sample in legal:
            logger.debug(
                "legal: (" + str(sample.replica) +
                "," + str(sample.trajectory) +
                "," + repr(sample.ensemble) +
                ")")
        # TODO: This can't go through numpy.random as Samples unwrap to
        # Snapshots when cast to arrays
        selected = random.choice(legal)
        logger.debug(
            "selected sample: (" + str(selected.replica) +
            "," + str(selected.trajectory) +
            "," + repr(selected.ensemble) +
            ")")
        return selected

    @abc.abstractmethod
    def move(self, sample_set):
        """
        Run the generation starting with the initial sample_set specified.

        Parameters
        ----------
        sample_set : SampleSet
            the initially used sampleset

        Returns
        -------
        samples : MoveChange
            the MoveChange instance describing the change from the old to
            the new SampleSet

        """

        return paths.EmptyMoveChange()  # pragma: no cover

    def __str__(self):
        if self.name == self.__class__.__name__:
            return self.__repr__()
        else:
            return self.name


class IdentityPathMover(PathMover):
    """
    The simplest Mover that does nothing !

    Notes
    -----
    Since is does nothing it is considered rejected everytime!
    It can be used to test function of PathMover

    Parameters
    ----------
    counts_as_trial : bool
        Whether this mover should count as a trial or not. If `True`, the
        `EmptyMoveChange` returned includes this mover, which means it gets
        counted as a trial in analysis of acceptance. If `False` (default),
        the mover for the returned move change is `None`, which does not get
        counted as a trial.
    """
    def __init__(self, counts_as_trial=False):
        super(IdentityPathMover, self).__init__()
        self.counts_as_trial = counts_as_trial

    def move(self, sample_set):
        mover = self if self.counts_as_trial else None
        return paths.EmptyMoveChange(mover=mover)


###############################################################################
# GENERATORS
###############################################################################

class SampleMover(PathMover):
    def __init__(self):
        super(SampleMover, self).__init__()

    def metropolis(self, trials):
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
        details : openpathsampling.Details
            Returns a Details object that contains information about the
            decision, i.e. total acceptance and random number

        """

        shoot_str = "MC in {cls} using samples {trials}"
        logger.info(shoot_str.format(cls=self.__class__.__name__,
                                     trials=trials))

        trial_dict = dict()
        for trial in trials:
            trial_dict[trial.ensemble] = trial

        accepted = True
        probability = 1.0

        # TODO: This isn't right. `bias` should be associated with the
        # change; not with each individual sample. ~~~DWHS
        for ens, sample in trial_dict.items():
            valid = ens(sample.trajectory, candidate=self._trust_candidate)
            if not valid:
                # one sample not valid reject
                accepted = False
                probability = 0.0
                break
            else:
                probability *= sample.bias

        rand = self._rng.random()

        if rand > probability:
            # rejected
            accepted = False

        details = {
            'metropolis_acceptance': probability,
            'metropolis_random': rand
        }

        if accepted:
            result_str = "accepted"
        else:
            result_str = ("rejected. Acceptance probabilty "
                          + str(probability))

        logger.info("Trial was " + result_str)

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

    def get_samples_from_sample_set(self, sample_set):
        """
        Select samples to use as input to the move core.

        See Also
        --------
        move_core
        move

        Parameters
        ----------
        sample_set : :class:`.SampleSet`
            current samples to use as potential input

        Returns
        -------
        list of :class:`.Sample`
            samples to use as input to the move core
        """
        ensembles = self._called_ensembles()
        samples = [self.select_sample(sample_set, ens) for ens in ensembles]
        return samples

    def move(self, sample_set):
        samples = self.get_samples_from_sample_set(sample_set)
        change = self.move_core(samples)
        return change

    def move_core(self, samples):
        """Core of the Monte Carlo move. Includes acceptance.

        See Also
        --------
        move

        Parameters
        ----------
        samples : list of :class:`.Sample`
            input samples from the correct ensembles of this object

        Returns
        -------
        :class:`.MoveChange`
            result MoveChange for this move
        """
        # this is separated out for reuse and remove dependence core MC move
        # dependence on the entire sample set (for parallelization)
        try:
            # pass these samples to the trial move which might throw
            # engine-specific exceptions if something goes wrong.
            # Most common should be `EngineNaNError` if nan is detected and
            # `EngineMaxLengthError`
            trials, call_details = self(*samples)

        except SampleNaNError as e:
            e.details.update({'rejection_reason': 'nan'})
            return paths.RejectedNaNSampleMoveChange(
                samples=e.trial_sample,
                mover=self,
                input_samples=samples,
                details=paths.Details(**e.details)
            )
        except SampleMaxLengthError as e:
            e.details.update({'rejection_reason': 'max_length'})
            return paths.RejectedMaxLengthSampleMoveChange(
                samples=e.trial_sample,
                mover=self,
                input_samples=samples,
                details=paths.Details(**e.details)
            )

        accepted, acceptance_details = self._accept(trials)

        # update details
        kwargs = {}
        kwargs.update(call_details)
        kwargs.update(acceptance_details)

        details = Details(**kwargs)

        # return change
        if accepted:
            return paths.AcceptedSampleMoveChange(
                samples=trials,
                mover=self,
                input_samples=samples,
                details=details
            )
        else:
            return paths.RejectedSampleMoveChange(
                samples=trials,
                mover=self,
                input_samples=samples,
                details=details
            )

    @abc.abstractmethod
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

    Notes
    -----

    A few comments for developers working with subclasses of
    ``EngineMover``: This class is intended to do most of the grunt work for
    a wide range of possible engine-based needs. Remember that your
    ``selector`` can select first or final points, e.g., to extend a move.
    In order to help you find your way through the ``EngineMover`` code,
    here is an overview of what various private methods do:

    * ``__call__``: Creates the trial. Two steps: (1) make the trajectory;
      (2) assemble a sample to return
    * ``_build_sample``: assembles the final sample
    * ``_make_forward_trajectory``/``_make_backward_trajectory``: creates
      the actual trajectory, using :class:`.PrefixTrajectoryEnsemble` or
      :class:`.SuffixTrajectoryEnsemble` to ensure reasonable behavior (see
      below for further discussion)
    * ``._run``: this is what is called by ``__call__``, and it in turn
      calls the functions to make the trajectories (depending on the nature
      of the mover). Frequently, this is the only thing to override (two-way
      shooting, shifting).
    """

    default_engine = None
    reject_max_length = True

    # this will store the engine attribute for all subclasses as well
    _included_attr = ['_engine']

    def __init__(self, ensemble, target_ensemble, selector, engine=None):
        super(EngineMover, self).__init__()
        self.selector = selector
        self.ensemble = ensemble
        self.target_ensemble = target_ensemble
        self._engine = engine
        self._trust_candidate = True  # can I safely do that?
        # I think that is safe. Note for future: if we come across a bug
        # based on this, an alternative would be to have the move strategy
        # set _trust_candidate when it builds the movers; that is likely to
        # be a little safer (although I think we can trust all candidates
        # from engine movers to actually be candidates)

    def to_dict(self):
        dct = super(EngineMover, self).to_dict()
        dct['engine'] = self.engine
        return dct

    @property
    def engine(self):
        if self._engine is not None:
            return self._engine
        else:
            return self.default_engine

    @engine.setter
    def engine(self, engine):
        self._engine = engine

    def _called_ensembles(self):
        return [self.ensemble]

    def _get_in_ensembles(self):
        return [self.ensemble]

    def _get_out_ensembles(self):
        return [self.target_ensemble]

    def __call__(self, input_sample):
        initial_trajectory = input_sample.trajectory
        shooting_index = self.selector.pick(initial_trajectory)

        try:
            trial_trajectory, run_details = self._run(initial_trajectory,
                                                      shooting_index)

        except paths.engines.EngineNaNError as e:
            trial, details = self._build_sample(
                input_sample, shooting_index, e.last_trajectory, 'nan')

            raise SampleNaNError('Sample with NaN', trial, details)

        except paths.engines.EngineMaxLengthError as e:
            trial, details = self._build_sample(
                input_sample, shooting_index, e.last_trajectory, 'max_length')

            if EngineMover.reject_max_length:
                raise SampleMaxLengthError('Sample with MaxLength', trial, details)

        else:
            trial, details = self._build_sample(
                input_sample, shooting_index, trial_trajectory)

        trials = [trial]
        details.update(run_details)

        return trials, details

    def _build_sample(
            self,
            input_sample,
            shooting_index,
            trial_trajectory,
            stopping_reason=None
            ):

        initial_trajectory = input_sample.trajectory

        if stopping_reason is None:
            bias = self.selector.probability_ratio(
                initial_trajectory[shooting_index],
                initial_trajectory,
                trial_trajectory
            )
        else:
            bias = 0.0

        # temporary test to make sure nothing went weird
        # old_bias = initial_point.sum_bias / trial_point.sum_bias
        # assert(abs(bias - old_bias) < 10e-6)
        # assert(initial_trajectory[shooting_index] in trial_trajectory)

        # we need to save the initial
        trial_details = {
            'initial_trajectory': initial_trajectory,
            'shooting_snapshot': initial_trajectory[shooting_index]
        }

        if stopping_reason is not None:
            trial_details['stopping_reason'] = stopping_reason

        trial = paths.Sample(
            replica=input_sample.replica,
            trajectory=trial_trajectory,
            ensemble=self.target_ensemble,
            parent=input_sample,
            mover=self,
            bias=bias
        )

        return trial, trial_details

    def _make_forward_trajectory(self, trajectory, shooting_index):
        initial_snapshot = trajectory[shooting_index]  # .copy()
        run_f = paths.PrefixTrajectoryEnsemble(self.target_ensemble,
                                               trajectory[0:shooting_index]
                                              ).can_append
        partial_trajectory = self.engine.generate(initial_snapshot,
                                                  running=[run_f])
        trial_trajectory = (trajectory[0:shooting_index] +
                            partial_trajectory)
        # TODO: this should check for overshoot; only works now if ensemble
        # doesn't overshoot
        return trial_trajectory

    def _make_backward_trajectory(self, trajectory, shooting_index):
        initial_snapshot = trajectory[shooting_index].reversed  # _copy()
        run_f = paths.SuffixTrajectoryEnsemble(self.target_ensemble,
                                               trajectory[shooting_index + 1:]
                                              ).can_prepend
        partial_trajectory = self.engine.generate(initial_snapshot,
                                                  running=[run_f])
        trial_trajectory = (partial_trajectory.reversed +
                            trajectory[shooting_index + 1:])
        # TODO: this should check for overshoot; only works now if ensemble
        # doesn't overshoot
        return trial_trajectory

    # direction is an abstract property to disallow instantiation
    # of the EngineMover unless we use a concrete subclass that sets this.
    # This is not super elegant but is the way to do it with abstract classes

    @abc.abstractproperty
    def direction(self):
        return 'unknown'

    def _run(self, trajectory, shooting_index):
        """Takes initial trajectory and shooting point; return trial
        trajectory"""
        shoot_str = "Running {sh_dir} from frame {fnum} in [0:{maxt}]"
        logger.info(shoot_str.format(
            fnum=shooting_index,
            maxt=len(trajectory) - 1,
            sh_dir=self.direction
        ))

        if self.direction == "forward":
            trial_trajectory = self._make_forward_trajectory(
                trajectory, shooting_index
            )
        elif self.direction == "backward":
            trial_trajectory = self._make_backward_trajectory(
                trajectory, shooting_index
            )
        else:
            raise RuntimeError("Unknown direction: " + str(self.direction))

        return trial_trajectory, {}


class ForwardShootMover(EngineMover):
    """A forward shooting sample generator
    """

    def __init__(self, ensemble, selector, engine=None):
        super(ForwardShootMover, self).__init__(
            ensemble=ensemble,
            target_ensemble=ensemble,
            selector=selector,
            engine=engine
        )

    @property
    def direction(self):
        return 'forward'


class BackwardShootMover(EngineMover):
    """A Backward shooting generator
    """

    def __init__(self, ensemble, selector, engine=None):
        super(BackwardShootMover, self).__init__(
            ensemble=ensemble,
            target_ensemble=ensemble,
            selector=selector,
            engine=engine
        )

    @property
    def direction(self):
        return 'backward'


class ForwardExtendMover(EngineMover):
    """
    A Sample Mover implementing Forward Extension
    """
    _direction = "forward"

    def __init__(self, ensemble, target_ensemble, engine=None):
        super(ForwardExtendMover, self).__init__(
            ensemble=ensemble,
            target_ensemble=target_ensemble,
            selector=paths.FinalFrameSelector(),
            engine=engine
        )

    @property
    def direction(self):
        return 'forward'


class BackwardExtendMover(EngineMover):
    """
    A Sample Mover implementing Backward Extension
    """
    _direction = "backward"

    def __init__(self, ensemble, target_ensemble, engine=None):
        super(BackwardExtendMover, self).__init__(
            ensemble=ensemble,
            target_ensemble=target_ensemble,
            selector=paths.FirstFrameSelector(),
            engine=engine
        )

    @property
    def direction(self):
        return 'backward'


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
        self._trust_candidate = True

        initialization_logging(logger=init_log, obj=self,
                               entries=['bias', 'ensemble1', 'ensemble2'])

    def _called_ensembles(self):
        return [self.ensemble1, self.ensemble2]

    def _get_in_ensembles(self):
        return [self.ensemble1, self.ensemble2]

    def _get_out_ensembles(self):
        return [self.ensemble1, self.ensemble2]

    def _generate_in_out(self):
        return InOutSet([
            InOut([
                ((self.ensemble1, self.ensemble2, 1 * InOut._use_move_type), 1),
                ((self.ensemble2, self.ensemble1, 1 * InOut._use_move_type), 1)
            ])
        ])

    def __call__(self, sample1, sample2):
        # convert sample to the language used here before
        trajectory1 = sample1.trajectory
        trajectory2 = sample2.trajectory
        ensemble1 = sample1.ensemble
        ensemble2 = sample2.ensemble
        replica1 = sample1.replica
        replica2 = sample2.replica

        from1to2 = ensemble2(trajectory1, candidate=self._trust_candidate)
        logger.debug("trajectory " + repr(trajectory1) +
                     " into ensemble " + repr(ensemble2) +
                     " : " + str(from1to2))
        from2to1 = ensemble1(trajectory2, candidate=self._trust_candidate)
        logger.debug("trajectory " + repr(trajectory2) +
                     " into ensemble " + repr(ensemble1) +
                     " : " + str(from2to1))

        trial1 = paths.Sample(
            replica=replica1,
            trajectory=trajectory1,
            ensemble=ensemble2,
            parent=sample1,
            mover=self
        )
        trial2 = paths.Sample(
            replica=replica2,
            trajectory=trajectory2,
            ensemble=ensemble1,
            parent=sample2,
            mover=self
        )

        return [trial1, trial2], {}


class StateSwapMover(SampleMover):
    def __init__(self, ensemble1, ensemble2, bias=None):
        """
        A move to swap states for state changing samples

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

    def _generate_in_out(self):
        return InOutSet([
            InOut([
                ((self.ensemble1, self.ensemble2, -1 * InOut._use_move_type), 1),
                ((self.ensemble2, self.ensemble1, -1 * InOut._use_move_type), 1)
            ])
        ])

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
            mover=self
        )
        trial2 = paths.Sample(
            replica=replica2,
            trajectory=trajectory2,
            ensemble=ensemble1,
            parent=sample2,
            mover=self
        )

        return [trial1, trial2], {}


###############################################################################
# SUBTRAJECTORY GENERATORS
###############################################################################

class SubtrajectorySelectMover(SampleMover):
    """
    Picks a subtrajectory satisfying the given subensemble.

    If there are no subtrajectories which satisfy the subensemble, this
    returns the zero-length trajectory.

    Attributes
    ----------
    ensemble : openpathsampling.Ensemble
        the set of allows samples to chose from
    sub_ensemble : openpathsampling.Ensemble
        the subensemble to be searched for
    n_l : int or None
        the number of subtrajectories that need to be found. If
        `None` every number of subtrajectories > 0 is okay.
        Otherwise the move is only accepted if exactly n_l subtrajectories
        are found.

    """

    _is_ensemble_change_mover = True

    def __init__(self, ensemble, sub_ensemble, n_l=None):
        super(SubtrajectorySelectMover, self).__init__(
        )
        self.n_l = n_l
        self.ensemble = ensemble
        self.sub_ensemble = sub_ensemble

    def _called_ensembles(self):
        return [self.ensemble]

    def _get_in_ensembles(self):
        return [self.ensemble]

    def _get_out_ensembles(self):
        return [self.sub_ensemble]

    @abc.abstractmethod
    def _choose(self, trajectory_list):
        return [], {}

    def __call__(self, trial):
        initial_trajectory = trial.trajectory
        replica = trial.replica
        logger.debug(
            "Working with replica " + str(replica) +
            " (" + str(initial_trajectory) + ")")

        subtrajs = self.sub_ensemble.split(initial_trajectory)
        logger.debug("Found " + str(len(subtrajs)) + " subtrajectories.")

        if (self.n_l is None and len(subtrajs) > 0) or \
                (self.n_l is not None and len(subtrajs) == self.n_l):
            subtraj, selection_details = self._choose(subtrajs)

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
            selection_details = {}

        return trials, selection_details


class RandomSubtrajectorySelectMover(SubtrajectorySelectMover):
    """
    Samples a random subtrajectory satisfying the given subensemble.

    If there are no subtrajectories which satisfy the subensemble, this
    returns the zero-length trajectory.

    Attributes
    ----------
    ensemble : openpathsampling.Ensemble
        the set of allows samples to chose from
    sub_ensemble : openpathsampling.Ensemble
        the subensemble to be searched for
    n_l : int or None
        the number of subtrajectories that need to be found. If
        `None` every number of subtrajectories > 0 is okay.
        Otherwise the move is only accepted if exactly n_l subtrajectories
        are found.

    """

    def _choose(self, trajectory_list):
        # Needed to prevent ragged nested sequence DeprWarning in np 1.20
        idx = self._rng.choice(len(trajectory_list))
        return trajectory_list[idx], {}


class FirstSubtrajectorySelectMover(SubtrajectorySelectMover):
    """
    Samples the first subtrajectory satifying the given subensemble.

    If there are no subtrajectories which satisfy the ensemble, this returns
    the zero-length trajectory.

    Attributes
    ----------
    ensemble : openpathsampling.Ensemble
        the set of allows samples to chose from
    sub_ensemble : openpathsampling.Ensemble
        the subensemble to be searched for
    n_l : int or None
        the number of subtrajectories that need to be found. If
        `None` every number of subtrajectories > 0 is okay.
        Otherwise the move is only accepted if exactly n_l subtrajectories
        are found.

    """

    def _choose(self, trajectory_list):
        return trajectory_list[0], {}


class FinalSubtrajectorySelectMover(SubtrajectorySelectMover):
    """
    Samples the final subtrajectory satifying the given subensemble.

    If there are no subtrajectories which satisfy the ensemble, this returns
    the zero-length trajectory.

    Attributes
    ----------
    ensemble : openpathsampling.Ensemble
        the set of allows samples to chose from
    sub_ensemble : openpathsampling.Ensemble
        the subensemble to be searched for
    n_l : int or None
        the number of subtrajectories that need to be found. If
        `None` every number of subtrajectories > 0 is okay.
        Otherwise the move is only accepted if exactly n_l subtrajectories
        are found.

    """

    def _choose(self, trajectory_list):
        return trajectory_list[-1], {}


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
        self._trust_candidate = True

    def _called_ensembles(self):
        return [self.ensemble]

    def _get_in_ensembles(self):
        return [self.ensemble]

    def _generate_in_out(self):
        return InOutSet([
            InOut([
                ((self.ensemble, self.ensemble, -1 * InOut._use_move_type), 1)
            ])
        ])

    def __call__(self, trial):
        trajectory = trial.trajectory
        ensemble = trial.ensemble
        replica = trial.replica

        reversed_trajectory = trajectory.reversed

        valid = ensemble(reversed_trajectory,
                         candidate=self._trust_candidate)
        logger.info("PathReversal move accepted: " + str(valid))

        bias = 1.0

        trial = paths.Sample(
            replica=replica,
            trajectory=reversed_trajectory,
            ensemble=ensemble,
            mover=self,
            parent=trial,
            bias=bias
        )

        return [trial], {}


class EnsembleHopMover(SampleMover):
    _is_ensemble_change_mover = True

    def __init__(
            self, ensemble, target_ensemble, change_replica=None, bias=None):
        """
        A Mover that allows the change between ensembles.

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

            .. code-block:: python

                {
                    'ensembles' : [ens_1, ens_2, ens_n],
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
        return [self.ensemble]

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

        logger.info(
            "Attempting ensemble hop from {e1} to {e2} replica ID {rid}".format(
                e1=repr(ens_from), e2=repr(ens_to), rid=repr(replica)))

        trajectory = rep_sample.trajectory
        logger.debug("  selected replica: " + str(replica))
        logger.debug("  initial ensemble: " + repr(rep_sample.ensemble))

        logger.info("Hop starts from legal ensemble: "+str(ens_from(trajectory)))
        logger.info("Hop ends in legal ensemble: "+str(ens_to(trajectory)))


        # TODO: remove this and generalize!!!
        if type(self.bias) is float:
            bias = self.bias
            logger.info("Using fixed bias " + str(bias))
        elif type(self.bias) is dict:
            # special dict
            ens = self.bias['ensembles']
            e1 = ens.index(ens_from)
            e2 = ens.index(ens_to)
            bias = float(self.bias['values'][e1, e2])
            logger.info("Using dict bias " + str(bias))
        else:
            bias = 1.0
            logger.info("Using default bias: self.bias == " + str(self.bias))

        trial = paths.Sample(
            replica=replica,
            trajectory=trajectory,
            ensemble=ens_to,
            mover=self,
            parent=rep_sample,
            bias=bias
        )

        details = {
            'initial_ensemble': ens_from,
            'trial_ensemble': ens_to,
            'bias': bias
        }

        return [trial], details


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

    def _generate_in_out(self):
        return InOutSet(set.union(*[sub.in_out for sub in self.submovers]))

    def _get_in_ensembles(self):
        return [sub.input_ensembles for sub in self.submovers]

    def _get_out_ensembles(self):
        return [sub.output_ensembles for sub in self.submovers]

    @abc.abstractmethod
    def _selector(self, sample_set):
        pass

    def select_mover(self, weights):
        p = np.array(weights)
        p /= sum(weights)
        logger.debug(self.name + " " + str(weights))
        idx = self._rng.choice(len(self.movers), p=p)

        logger_str = "{name} ({cls}) selecting {mtype} (index {idx})"
        logger.info(logger_str.format(
            name=self.name,
            cls=self.__class__.__name__,
            idx=idx,
            mtype=self.movers[idx].name
        ))

        mover = self.movers[idx]

        kwargs = {
            'choice': idx,
            'chosen_mover': mover,
            'probability': weights[idx] / sum(weights),
            'weights': weights
        }

        details = Details(**kwargs)
        return mover, details

    def move(self, sample_set):
        weights = self._selector(sample_set)
        mover, details = self.select_mover(weights)
        subchange = mover.move(sample_set)

        path = paths.RandomChoiceMoveChange(
            subchange=subchange,
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

    def _selector(self, sample_set):
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

    def _selector(self, sample_set):
        if self.weights is None:
            weights = [1.0] * len(self.movers)
        else:
            weights = list(self.weights)  # make a copy

        # this is implemented by setting all weights locally to zero that
        # correspond to movers that will potentially fail since the required
        # input ensembles are not present in the sample_set

        present_ensembles = sample_set.ensembles

        for idx, mover in enumerate(self.movers):
            for ens in mover.input_ensembles:
                if ens not in present_ensembles:
                    # ens might be required but is not present
                    weights[idx] = 0.0

        return weights


class FirstAllowedMover(SelectionMover):
    """
    Chooses a first mover that has samples in all required ensembles.

    A mover can only safely be run, if all inputs can be satisfied.
    This will pick the first mover from the list where all ensembles
    from input_ensembles are found.

    Attributes
    ----------
    movers : list of PathMover
        the PathMovers to choose from
    weights : list of floats
        the relative weight of each PathMover (does not need to be normalized)

    """

    def _selector(self, sample_set):
        weights = [1.0] * len(self.movers)

        present_ensembles = sample_set.ensembles

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

    A mover can only safely be run, if all inputs can be satisfied.
    This will pick the last mover from the list where all ensembles
    from input_ensembles are found.

    Attributes
    ----------
    movers : list of PathMover
        the PathMovers to choose from
    weights : list of floats
        the relative weight of each PathMover (does not need to be normalized)

    """

    def _selector(self, sample_set):
        weights = [1.0] * len(self.movers)

        present_ensembles = sample_set.ensembles

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

    Returns a SequentialMoveChange of the if_move movepath and the then_move
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

    def _generate_in_out(self):
        return InOutSet({
            self.if_mover.in_out + self.then_mover.in_out
        } | {
            self.if_mover.in_out + self.else_mover.in_out
        })

    def sub_replica_state(self, replica_states):
        if_replica_states = self.if_mover.in_out.move(replica_states)
        return [
            if_replica_states,
            self.then_mover.in_out.move(if_replica_states),
            self.else_mover.in_out.move(if_replica_states),
        ]

    @property
    def submovers(self):
        return [self.if_mover, self.then_mover, self.else_mover]

    def _get_in_ensembles(self):
        return [sub.input_ensembles for sub in self.submovers]

    def _get_out_ensembles(self):
        return [sub.output_ensembles for sub in self.submovers]

    def move(self, sample_set):
        subglobal = sample_set

        ifclause = self.if_mover.move(subglobal)
        samples = ifclause.results
        subglobal = subglobal.apply_samples(samples)

        if ifclause.accepted:
            if self.then_mover is not None:
                resultclause = self.then_mover.move(subglobal)
            else:
                resultclause = paths.EmptyMoveChange()
        else:
            if self.else_mover is not None:
                resultclause = self.else_mover.move(subglobal)
            else:
                resultclause = paths.EmptyMoveChange()

        return paths.SequentialMoveChange([ifclause, resultclause], mover=self)


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

    def _generate_in_out(self):
        total = self.submovers[0].in_out

        for pp in range(1, len(self.submovers)):
            new_pos = self.submovers[pp].in_out
            total = InOutSet(set.union(*[
                total,
                InOutSet(
                    sum([total, new_pos], InOutSet())
                ),
                new_pos
            ]))

        return total

    def sub_replica_state(self, replica_states):
        ret = list()
        ret.append(replica_states)
        for sub in self.submovers[:-1]:
            replica_states = sub.in_out.move(replica_states)
            ret.append(replica_states)

        return ret

    def _get_in_ensembles(self):
        return [sub.input_ensembles for sub in self.submovers]

    def _get_out_ensembles(self):
        return [sub.output_ensembles for sub in self.submovers]

    def move(self, sample_set):
        logger.debug("Starting sequential move")

        subglobal = sample_set
        movechanges = []

        for mover in self.movers:
            logger.debug("Starting sequential move step " + str(mover))

            # Run the sub mover
            movepath = mover.move(subglobal)
            samples = movepath.results
            subglobal = subglobal.apply_samples(samples)
            movechanges.append(movepath)

        return paths.SequentialMoveChange(movechanges, mover=self)


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

    def _generate_in_out(self):
        total = self.submovers[0].in_out

        for pp in range(1, len(self.submovers)):
            new_pos = self.submovers[pp].in_out
            total = InOutSet(set.union(*[
                total,
                InOutSet(
                    sum([total, new_pos], InOutSet())
                )
            ]))

        return total

    def move(self, sample_set):
        logger.debug("==== BEGINNING " + self.name + " ====")
        subglobal = paths.SampleSet(sample_set)
        movechanges = []
        for mover in self.movers:
            logger.info(str(self.name)
                        + " starting mover index " + str(self.movers.index(mover))
                        + " (" + mover.name + ")"
                       )
            # Run the sub mover
            movepath = mover.move(subglobal)
            samples = movepath.results
            subglobal = subglobal.apply_samples(samples)
            movechanges.append(movepath)
            if not movepath.accepted:
                break

        logger.debug("==== FINISHING " + self.name + " ====")
        return paths.PartialAcceptanceSequentialMoveChange(
            movechanges, mover=self)


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

    def _generate_in_out(self):
        return InOutSet(sum([sub.in_out for sub in self.submovers], InOutSet()))

    def move(self, sample_set):
        logger.debug("Starting conditional sequential move")

        subglobal = sample_set
        movechanges = []

        for mover in self.movers:
            logger.debug("Starting sequential move step " + str(mover))

            # Run the sub mover
            movepath = mover.move(subglobal)
            samples = movepath.results
            subglobal = subglobal.apply_samples(samples)
            movechanges.append(movepath)

            if not movepath.accepted:
                break

        return paths.ConditionalSequentialMoveChange(
            movechanges, mover=self)


class NonCanonicalConditionalSequentialMover(ConditionalSequentialMover):
    """ Special mover for reactive flux simulation.

    This mover inherits from :class:`.ConditionalSequentialMover` and
    alters only the `move` method to return the output of the corresponding
    :class:`.NonCanonicalConditionalSequentialMoveChange`.
    """
    _is_canonical = False

    def move(self, sample_set):
        change = super(NonCanonicalConditionalSequentialMover,
                       self).move(sample_set)
        return paths.NonCanonicalConditionalSequentialMoveChange(
            subchanges=change.subchanges,
            mover=change.mover,
            details=change.details
        )


# class ReplicaIDChangeMover(PathMover):
#     """
#     Changes the replica ID for a path.
#     """
#
#     def __init__(self, replica_pair):
#         super(ReplicaIDChangeMover, self).__init__()
#         self.replica_pair = replica_pair
#         initialization_logging(logger=init_log, obj=self,
#                                entries=['replica_pairs'])
#
#     def move(self, sample_set):
#         rep_from = self.replica_pair[0]
#         rep_to = self.replica_pair[1]
#         rep_sample = self.select_sample(sample_set,
#                                         replicas=rep_from)
#
#         logger.info(
#             "Creating new sample from replica ID " + str(rep_from) +
#             " and putting it in replica ID " + str(rep_to))
#
#         # note: currently this clones into a new replica ID. We might later
#         # want to kill the old replica ID (and possibly rename this mover).
#
#         new_sample = paths.Sample(
#             replica=rep_to,
#             ensemble=rep_sample.ensemble,
#             trajectory=rep_sample.trajectory,
#             parent=rep_sample,
#             mover=self
#         )
#
#         # Can be used to remove the old sample. Not used yet!
#         # kill_sample = paths.Sample(
#         #     replica=rep_from,
#         #     trajectory=None,
#         #     ensemble=rep_sample.ensemble,
#         #     parent=None,
#         #     mover=self
#         # )
#
#         kwargs = {
#             'rep_from': rep_from,
#             'rep_to': rep_to
#         }
#
#         details = Details(**kwargs)
#
#         return paths.AcceptedSampleMoveChange(
#             samples=[new_sample],
#             mover=self,
#             input_samples=[rep_sample],
#             details=details
#         )


class SubPathMover(PathMover):
    """Mover that delegates to a single submover
    """

    def __init__(self, mover):
        """
        Parameters
        ----------
        mover : :class:`openpathsampling.PathMover`
            the submover to be delegated to
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

    def _generate_in_out(self):
        return self.mover.in_out

    def sub_replica_state(self, replica_states):
        return [replica_states]

    def move(self, sample_set):
        subchange = self.mover.move(sample_set)
        change = paths.SubMoveChange(
            subchange=subchange,
            mover=self
        )
        return change


class EnsembleFilterMover(SubPathMover):
    """Mover that return only samples from specified ensembles
    """

    def __init__(self, mover, ensembles):
        """
        Parameters
        ----------
        mover : :class:`openpathsampling.PathMover`
            the submover to be delegated to
        ensembles : nested list of :class:`openpathsampling.Ensemble` or None
            the ensemble specification
        """
        super(EnsembleFilterMover, self).__init__(mover)
        self.ensembles = ensembles

        if not set(self.mover.output_ensembles) & set(self.ensembles):
            # little sanity check, if the underlying move will be removed by the
            # filter throw a warning
            raise ValueError(
                'Your filter removes the underlying move completely. ' +
                'Please check your ensembles and submovers!')

    def move(self, sample_set):
        # TODO: This will only pass filtered samples. We might split
        # this into an separate input and output filter if only one
        # side is needed

        filtered_globalstate = paths.SampleSet([
            samp for samp in sample_set if samp.ensemble in self.ensembles
        ])
        subchange = self.mover.move(filtered_globalstate)
        change = paths.FilterByEnsembleMoveChange(
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

    def _generate_in_out(self):
        return self.mover.in_out.filter(self.ensembles)

    def sub_replica_state(self, replica_states):
        return [{rs.filter(self.ensembles) for rs in replica_states}]


class SpecializedRandomChoiceMover(RandomChoiceMover):
    """
    Superclass for movers that are random choice between two SampleMovers

    This requires that all submovers accept the same list of samples.
    """
    @classmethod
    def from_dict(cls, dct):
        mover = cls.__new__(cls)
        super(cls, mover).__init__(movers=dct['movers'])
        return mover

    def to_dict(self):
        dct = super(SpecializedRandomChoiceMover, self).to_dict()
        dct['movers'] = self.movers
        return dct

    def move_core(self, samples):
        weights = self.weights
        mover, details = self.select_mover(weights)
        subchange = mover.move_core(samples)
        change = paths.RandomChoiceMoveChange(
            subchange=subchange,
            mover=self,
            details=details
        )
        return change

class OneWayShootingMover(SpecializedRandomChoiceMover):
    """
    One-way (stochastic) shooting mover

    OneWayShootingMover is a special case of a RandomChoiceMover which
    combines gives a 50/50 chance of selecting either a ForwardShootMover or
    a BackwardShootMover. Both submovers use the same shooting point
    selector, and both apply to the same ensembles and replicas.

    Attributes
    ----------
    selector : :class:`openpathsampling.ShootingPointSelector`
        The shooting point selection scheme
    ensemble : :class:`openpathsampling.Ensemble`
        Ensemble for this shooting mover
    """

    def __init__(self, ensemble, selector, engine=None):
        movers = [
            ForwardShootMover(ensemble=ensemble,
                              selector=selector,
                              engine=engine),
            BackwardShootMover(ensemble=ensemble,
                               selector=selector,
                               engine=engine)
        ]
        super(OneWayShootingMover, self).__init__(movers=movers)

    @classmethod
    def from_dict(cls, dct):
        mover = cls.__new__(cls)
        # override with stored movers and use the init of the super class
        # this assumes that the super class has movers as its signature
        super(cls, mover).__init__(movers=dct['movers'])
        return mover

    @property
    def ensemble(self):
        return self.movers[0].ensemble

    @property
    def selector(self):
        return self.movers[0].selector

    @property
    def engine(self):
        return self.movers[0].engine


class OneWayExtendMover(SpecializedRandomChoiceMover):
    """
    OneWayShootingMover is a special case of a RandomChoiceMover which
    gives a 50/50 chance of selecting either a ForwardExtendMover or
    a BackwardExtendMover. Both submovers use the same same ensembles
    and replicas.

    Attributes
    ----------
    ensemble : :class:`openpathsampling.Ensemble`
        valid ensemble
    """

    def __init__(self, ensemble, target_ensemble, engine=None):
        movers = [
            ForwardExtendMover(ensemble=ensemble,
                               target_ensemble=target_ensemble,
                               engine=engine),
            BackwardExtendMover(ensemble=ensemble,
                                target_ensemble=target_ensemble,
                                engine=engine)
        ]
        super(OneWayExtendMover, self).__init__(movers=movers)

    @property
    def engine(self):
        return self.movers[0].engine


class AbstractTwoWayShootingMover(EngineMover):
    def __init__(self, ensemble, selector, modifier, engine=None):
        super(AbstractTwoWayShootingMover, self).__init__(
            ensemble=ensemble,
            target_ensemble=ensemble,
            selector=selector,
            engine=engine
        )
        self.modifier = modifier

    # required for concrete class; not really used
    @property
    def direction(self):  # pragma: no cover
        return 'bidrectional'

    def _make_forward_trajectory(self, trajectory, initial_snapshot,
                                 shooting_index):
        fwd_ens = paths.PrefixTrajectoryEnsemble(
            self.target_ensemble,
            trajectory[0:shooting_index]
        )
        fwd_partial = self.engine.generate(initial_snapshot,
                                           running=[fwd_ens.can_append])
        return fwd_partial

    def _make_backward_trajectory(self, trajectory, initial_snapshot,
                                  shooting_index):
        # run backward
        bkwd_ens = paths.SuffixTrajectoryEnsemble(
            self.target_ensemble,
            trajectory[shooting_index + 1:]
        )
        bkwd_partial = self.engine.generate(initial_snapshot.reversed,
                                            running=[bkwd_ens.can_prepend])
        return bkwd_partial

    def _run(self, trajectory, shooting_index):
        # to override the default implementation in EngineMover
        raise NotImplementedError

class ForwardFirstTwoWayShootingMover(AbstractTwoWayShootingMover):
    def _run(self, trajectory, shooting_index):
        """
        The actual shooting process (after shooting point is chosen).

        Parameters
        ----------
        trajectory : :class:`.Trajectory`
            input trajectory
        shooting_index : int
            index of the shooting point within `trajectory`

        Returns
        -------
        trial_trajectory : :class:`.Trajectory`
            the resulting trial trajectory
        details : dict
            details dictionary (includes modified shooting point)
        """
        shoot_str = "Running {sh_dir} from frame {fnum} in [0:{maxt}]"
        logger.info(shoot_str.format(
            fnum=shooting_index,
            maxt=len(trajectory) - 1,
            sh_dir="Forward-first"
        ))

        original = trajectory[shooting_index]
        modified = self.modifier(original)

        fwd_partial = self._make_forward_trajectory(trajectory, modified,
                                                    shooting_index)
        # TODO: come up with a test that shows why you need mid_traj here;
        # should be a SeqEns with OptionalEnsembles. Exact example is hard!
        mid_traj = trajectory[0:shooting_index] + fwd_partial
        bkwd_partial = self._make_backward_trajectory(mid_traj, modified,
                                                      shooting_index)

        # join the two
        trial_trajectory = bkwd_partial.reversed + fwd_partial[1:]

        details = {'modified_shooting_snapshot': modified}
        return trial_trajectory, details


class BackwardFirstTwoWayShootingMover(AbstractTwoWayShootingMover):
    def _run(self, trajectory, shooting_index):
        """
        The actual shooting process (after shooting point is chosen).

        Parameters
        ----------
        trajectory : :class:`.Trajectory`
            input trajectory
        shooting_index : int
            index of the shooting point within `trajectory`

        Returns
        -------
        trial_trajectory : :class:`.Trajectory`
            the resulting trial trajectory
        details : dict
            details dictionary (includes modified shooting point)
        """
        shoot_str = "Running {sh_dir} from frame {fnum} in [0:{maxt}]"
        logger.info(shoot_str.format(
            fnum=shooting_index,
            maxt=len(trajectory) - 1,
            sh_dir="Backward-first"
        ))

        original = trajectory[shooting_index]
        modified = self.modifier(original)

        bkwd_partial = self._make_backward_trajectory(trajectory, modified,
                                                      shooting_index)
        #logger.info("Complete backward shot (length " +
                    #str(len(bkwd_partial)) + ")")
        # TODO: come up with a test that shows why you need mid_traj here;
        # should be a SeqEns with OptionalEnsembles. Exact example is hard!
        mid_traj = bkwd_partial.reversed + trajectory[shooting_index + 1:]
        mid_traj_shoot_idx = len(bkwd_partial) - 1
        fwd_partial = self._make_forward_trajectory(mid_traj, modified,
                                                    mid_traj_shoot_idx)
        #logger.info("Complete forward shot (length " +
                    #str(len(fwd_partial)) + ")")

        # join the two
        trial_trajectory = bkwd_partial.reversed + fwd_partial[1:]

        details = {'modified_shooting_snapshot': modified}
        return trial_trajectory, details


class TwoWayShootingMover(SpecializedRandomChoiceMover):
    def __init__(self, ensemble, selector, modifier, engine=None):
        movers = [
            ForwardFirstTwoWayShootingMover(
                ensemble=ensemble,
                selector=selector,
                modifier=modifier,
                engine=engine
            ),
            BackwardFirstTwoWayShootingMover(
                ensemble=ensemble,
                selector=selector,
                modifier=modifier,
                engine=engine
            )
        ]
        super(TwoWayShootingMover, self).__init__(movers=movers)

    @property
    def ensemble(self):
        return self.movers[0].ensemble

    @property
    def selector(self):
        return self.movers[0].selector

    @property
    def modifier(self):
        return self.movers[0].modifier


class MinusMover(SubPathMover):
    """
    Instance of a MinusMover.

    The minus move combines a replica exchange with path extension to swap
    paths between the innermost regular TIS interface ensemble and the minus
    interface ensemble. This is particularly useful for improving sampling
    of path space.
    """
    _is_canonical = True

    def __init__(self, minus_ensemble, innermost_ensembles, engine=None):
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
        sub_trajectory_selector.named("MinusSubtrajectoryChooser")

        repexs = [ReplicaExchangeMover(
            ensemble1=segment,
            ensemble2=inner
        ) for inner in innermost_ensembles]

        repex_chooser = RandomChoiceMover(repexs)
        repex_chooser.named("InterfaceSetChooser")

        extension_mover = RandomChoiceMover([
            ForwardExtendMover(
                ensemble=segment,
                target_ensemble=minus_ensemble,
                engine=engine
            ),
            BackwardExtendMover(
                ensemble=segment,
                target_ensemble=minus_ensemble,
                engine=engine
            )
        ])

        extension_mover.named("MinusExtensionDirectionChooser")
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

    def move(self, sample_set):
        change = super(MinusMover, self).move(sample_set)
        cond_seq_changes = change.subchanges[0].subchanges[0].subchanges
        seg_swap = None
        if len(cond_seq_changes) >= 2:
            seg_swap = cond_seq_changes[1].subchanges[0].trials

        ext_traj = None
        if len(cond_seq_changes) >= 3:
            ext_traj = cond_seq_changes[2].subchanges[0].trials[0].trajectory

        details = Details(segment_swap_samples=seg_swap,
                          extension_trajectory=ext_traj)
        if change.details is None:
            change.details = details

        return change

class SingleReplicaMinusMover(MinusMover):
    """
    Minus mover for single replica TIS.

    In SRTIS, the minus mover doesn't actually keep an active sample in the
    minus interface. Instead, it just puts the newly generated segment into
    the innermost ensemble.
    """

    def __init__(self, minus_ensemble, innermost_ensembles,
                 bias=None, engine=None):
        try:
            innermost_ensembles = list(innermost_ensembles)
        except TypeError:
            innermost_ensembles = [innermost_ensembles]

        if bias is None:
            bias = ""  # TODO temp for storage until real bias
        self.bias = bias
        self.minus_ensemble = minus_ensemble
        self.innermost_ensembles = innermost_ensembles

        # TODO: Until we have automated detailed balance calculations, I
        # think this will only be valid in the case of only one innermost
        # ensemble.  But I think you only want to use it in the case of only
        # one innermost ensemble anyway. The following warns us:
        if len(innermost_ensembles) > 1:
            logger.warning(
                "Probably shouldn't use SingleReplicaMinusMover with MISTIS")

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
            ForwardExtendMover(segment, minus_ensemble, engine=engine),
            FinalSubtrajectorySelectMover(minus_ensemble, segment),
            hop_segment_to_innermost
        ])

        backward_minus = ConditionalSequentialMover([
            hop_innermost_to_segment,
            BackwardExtendMover(segment, minus_ensemble, engine=engine),
            FirstSubtrajectorySelectMover(minus_ensemble, segment),
            hop_segment_to_innermost
        ])

        mover = EnsembleFilterMover(RandomChoiceMover([backward_minus,
                                                       forward_minus]),
                                    ensembles=innermost_ensembles)

        # we skip MinusMover's init and go to the grandparent
        super(MinusMover, self).__init__(mover)

    def move(self, sample_set):
        # skip the MinusMover's implementation
        return super(MinusMover, self).move(sample_set)


class PathSimulatorMover(SubPathMover):
    """
    This just wraps a mover and references the used pathsimulator
    """

    def __init__(self, mover, pathsimulator):
        super(PathSimulatorMover, self).__init__(mover)
        self.pathsimulator = pathsimulator

    def move(self, sample_set, step=-1):
        details = Details(
            step=step
        )

        return paths.PathSimulatorMoveChange(
            self.mover.move(sample_set),
            mover=self,
            details=details
        )


def PathReversalSet(ensembles):
    return list(map(PathReversalMover, ensembles))


class Details(StorableObject):
    """Details of an object. Can contain any data
    """

    def __init__(self, **kwargs):
        super(Details, self).__init__()
        for key, value in kwargs.items():
            setattr(self, key, value)

    _print_repr_types = [paths.Ensemble]
    _print_nothing_keys = ["__uuid__"]

    def __str__(self):
        # primarily for debugging/interactive use
        mystr = ""
        for key in self.__dict__.keys():
            obj = self.__dict__[key]
            if key in self._print_nothing_keys:
                pass  # do nothing!
            elif any([isinstance(obj, tt) for tt in self._print_repr_types]):
                mystr += str(key) + " = " + repr(obj) + '\n'
            else:
                mystr += str(key) + " = " + str(self.__dict__[key]) + '\n'
        return mystr


@has_deprecations
@deprecate(MOVE_DETAILS)
class MoveDetails(Details):
    """Details of the move as applied to a given replica

    Specific move types may have add several other attributes for each
    MoveDetails object. For example, shooting moves will also include
    information about the shooting point selection, etc.
    """

    def __init__(self, **kwargs):
        super(MoveDetails, self).__init__(**kwargs)


# leave this for potential backwards compatibility
@has_deprecations
@deprecate(SAMPLE_DETAILS)
class SampleDetails(Details):
    """Details of a sample
    """

    def __init__(self, **kwargs):
        super(SampleDetails, self).__init__(**kwargs)
