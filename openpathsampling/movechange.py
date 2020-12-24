__author__ = 'Jan-Hendrik Prinz'

import logging

import openpathsampling as paths
from openpathsampling.netcdfplus import StorableObject, lazy_loading_attributes
from openpathsampling.netcdfplus import DelayedLoader
from .treelogic import TreeMixin

logger = logging.getLogger(__name__)


# @lazy_loading_attributes('details')
class MoveChange(TreeMixin, StorableObject):
    '''
    A class that described the concrete realization of a PathMove.

    Attributes
    ----------
    mover : PathMover
        The mover that generated this MoveChange
    samples : list of Sample
        A list of newly generated samples by this particular move.
        Only used by node movers like RepEx or Shooters
    subchanges : list of MoveChanges
        the MoveChanges created by submovers
    details : Details
        an object that contains MoveType specific attributes and information.
        E.g. for a RandomChoiceMover which Mover was selected.
    '''

    details = DelayedLoader()

    def __init__(self, subchanges=None, samples=None, mover=None,
                 details=None, input_samples=None):
        StorableObject.__init__(self)

        self._lazy = {}

        self._len = None
        self._collapsed = None
        self._results = None
        self._trials = None
        self._accepted = None
        self.mover = mover
        if subchanges is None:
            self.subchanges = []
        else:
            self.subchanges = subchanges

        if samples is None:
            self.samples = []
        else:
            self.samples = samples

        if input_samples is None:
            self.input_samples = []
        else:
            self.input_samples = input_samples

        self.details = details

    def __getattr__(self, item):
        # try to get attributes from details dict
        try:
            return getattr(self.details, item)
        except AttributeError as e:
            msg = "{0} not found in change's details".format(str(item))
            if not e.args:
                e.args = [msg]
            else:
                e.args = tuple([e.args[0] + "; " + msg] + list(e.args[1:]))
            raise

    def to_dict(self):
        return {
            'mover': self.mover,
            'details': self.details,
            'samples': self.samples,
            'input_samples': self.input_samples,
            'subchanges': self.subchanges,
            'cls': self.__class__.__name__
        }

    # hook for TreeMixin
    @property
    def _subnodes(self):
        return self.subchanges

    @property
    def submovers(self):
        return [ch.mover for ch in self.subchanges]

    @property
    def subchange(self):
        """
        Return the single/only sub-movechange if there is only one.

        Returns
        -------
        MoveChange
        """
        if len(self.subchanges) == 1:
            return self.subchanges[0]
        else:
            # TODO: might raise exception
            return None

    @staticmethod
    def _default_match(original, test):
        if isinstance(test, paths.MoveChange):
            return original is test
        elif isinstance(test, paths.PathMover):
            return original.mover is test
        elif issubclass(test, paths.PathMover):
            return original.mover.__class__ is test
        else:
            return False

    def movetree(self):
        """
        Return a tree with the movers of each node

        Notes
        -----
        This is equivalent to
        `tree.map_tree(lambda x : x.mover)`
        """
        return self.map_tree(lambda x: x.mover)

    @property
    def identifier(self):
        return self.mover

    @property
    def collapsed_samples(self):
        """
        Return a collapsed set of samples with non used samples removed

        This is the minimum required set of samples to keep the `MoveChange`
        correct and allow to target sampleset to be correctly created.
        These are the samples used by `.closed`

        Examples
        --------
        Assume that you run 3 shooting moves for replica #1. Then only the
        last of the three really matters for the target sample_set since #1
        will be replaced by #2 which will be replaced by #3. So this function
        will return only the last sample.

        """
        if self._collapsed is None:
            s = paths.SampleSet([]).apply_samples(self.results)

            # keep order just for being thorough
            self._collapsed = [
                samp for samp in self.results
                if samp in s
            ]

        return self._collapsed

    @property
    def accepted(self):
        """
        Returns if this particular move was accepted.

        Mainly used for rejected samples.

        Notes
        -----
        Acceptance is determined from the number of resulting samples. If at
        least one sample is returned then this move will change the sampleset
        and is considered an accepted change.
        """
        if self._accepted is None:
            self._accepted = len(self.results) > 0

        return self._accepted

    def __add__(self, other):
        """
        This allows to use `+` to create SequentialPMCs

        Notes
        -----
        You can also use this to apply several changes
        >>> new_sset = old_sset + change1 + change2
        >>> new_sset = old_sset + (change1 + change2)
        """
        if isinstance(other, MoveChange):
            return SequentialMoveChange([self, other])
        else:
            raise ValueError('Only MoveChanges can be combined')

    @property
    def results(self):
        """
        Returns a list of all samples that are accepted in this move

        This contains unnecessary, but accepted samples, too.

        Returns
        -------
        list of Samples
            the list of samples that should be applied to the SampleSet
        """
        if self._results is None:
            self._results = self._get_results()

        return self._results

    def _get_results(self):
        """
        Determines all relevant accepted samples for this move

        Includes all accepted samples also from subchanges

        Returns
        -------
        list of Sample
            the list of accepted samples for this move
        """
        return []

    @property
    def trials(self):
        """
        Returns a list of all samples generated during the PathMove.

        This includes all accepted and rejected samples (which does NOT
        include hidden samples yet)

        """
        if self._trials is None:
            self._trials = self._get_trials()
        return self._trials

    def _get_trials(self):
        """
        Determines all samples for this move

        Includes all samples also from subchanges

        Returns
        -------
        list of Sample
            the list of all samples generated for this move

        Notes
        -----
        This function needs to be implemented for custom changes
        """
        return []

    def __str__(self):
        if self.accepted:
            return 'SampleMove : %s : %s : %d samples' % (self.mover.cls, self.accepted, len(self.trials)) + ' ' + str(self.trials) + ''
        else:
            return 'SampleMove : %s : %s :[]' % (self.mover.cls, self.accepted)

    @property
    def canonical(self):
        """
        Return the first non single-subchange

        Notes
        -----
        Usually a mover that returns a single subchange is for deciding what to
        do rather than describing what is actually happening. This property
        returns the first mover that is not one of these delegating movers and
        contains information of what has been done in this move.

        What you are usually interested in is `.canonical.mover` to get the
        relevant mover.

        Examples
        --------
        >>> a = OnewayShootingMover()
        >>> change = a.move(sset)
        >>> change.canonical.mover  # returns either Forward or Backward
        """
        pmc = self
        while pmc.subchange is not None:
            if pmc.mover.is_canonical is True:
                return pmc
            pmc = pmc.subchange

        return pmc

    @property
    def description(self):
        """
        Return a compact representation of the change
        """
        subs = self.subchanges
        if len(subs) == 0:
            return str(self.mover)
        elif len(subs) == 1:
            return subs[0].description
        else:
            return ':'.join([sub.description for sub in subs])


class EmptyMoveChange(MoveChange):
    """
    A MoveChange representing no changes
    """
    def __init__(self, mover=None, details=None):
        super(EmptyMoveChange, self).__init__(mover=mover, details=details)

    def __str__(self):
        return ''

    def _get_trials(self):
        return []

    def _get_results(self):
        return []


class SampleMoveChange(MoveChange):
    """
    A MoveChange representing the application of samples.

    This is the most common MoveChange and all other moves use this
    as leaves and on the lowest level consist only of `SampleMoveChange`
    """
    def __init__(self, samples, mover=None, details=None, input_samples=None):
        """
        Parameters
        ----------
        samples : list of Samples
            a list of trial samples that are used in this change
        mover : PathMover
            the generating PathMover
        details : Details
            a details object containing specifics about the change

        Attributes
        ----------
        samples
        mover
        details
        """
        super(SampleMoveChange, self).__init__(
            mover=mover,
            details=details,
            input_samples=input_samples
        )

        if samples.__class__ is paths.Sample:
            samples = [samples]

        self.samples = samples

    def _get_results(self):
        return []

    def _get_trials(self):
        return self.samples


class AcceptedSampleMoveChange(SampleMoveChange):
    """
    Represents an accepted SamplePMC

    This will return the trial samples also as its result, hence it is
    accepted.
    """
    def _get_trials(self):
        return self.samples

    def _get_results(self):
        return self.samples


class RejectedSampleMoveChange(SampleMoveChange):
    """
    Represents an rejected SamplePMC

    This will return no samples also as its result, hence it is
    rejected.
    """


class RejectedNaNSampleMoveChange(RejectedSampleMoveChange):
    """
    Represents an rejected SamplePMC because of occurance of NaN

    This will return no samples as its result, hence it is rejected.
    """
    pass


class RejectedMaxLengthSampleMoveChange(RejectedSampleMoveChange):
    """
    Represents an rejected SamplePMC because of hitting the max length limit

    This will return no samples as its result, hence it is rejected.
    """
    pass


class SequentialMoveChange(MoveChange):
    """
    SequentialMoveChange has no own samples, only inferred Sampled from the
    underlying MovePaths
    """
    def __init__(self, subchanges, mover=None, details=None):
        """
        Parameters
        ----------
        subchanges : list of MoveChanges
            a list of MoveChanges to be applied in sequence
        mover
        details

        Attributes
        ----------
        subchanges
        mover
        details
        """
        super(SequentialMoveChange, self).__init__(mover=mover, details=details)
        self.subchanges = subchanges

    def _get_results(self):
        samples = []
        for subchange in self.subchanges:
            samples = samples + subchange.results
        return samples

    def _get_trials(self):
        samples = []
        for subchange in self.subchanges:
            samples = samples + subchange.trials
        return samples

    def __str__(self):
        return 'SequentialMove : %s : %d samples\n' % \
               (self.accepted, len(self.results)) + \
               MoveChange._indent('\n'.join(map(str, self.subchanges)))


class PartialAcceptanceSequentialMoveChange(SequentialMoveChange):
    """
    PartialAcceptanceSequentialMovePath has no own samples, only inferred
    Sampled from the underlying MovePaths
    """

    def _get_results(self):
        changes = []
        for subchange in self.subchanges:
            if subchange.accepted:
                changes.extend(subchange.results)
            else:
                break

        return changes

    def __str__(self):
        return 'PartialAcceptanceMove : %s : %d samples\n' % \
               (self.accepted, len(self.results)) + \
               MoveChange._indent('\n'.join(map(str, self.subchanges)))


class ConditionalSequentialMoveChange(SequentialMoveChange):
    """
    ConditionalSequentialMovePath has no own samples, only inferred Samples
    from the underlying MovePaths
    """

    def _get_results(self):
        changes = []
        for subchange in self.subchanges:
            if subchange.accepted:
                changes.extend(subchange.results)
            else:
                return []

        return changes

    def __str__(self):
        return 'ConditionalSequentialMove : %s : %d samples\n' % \
               (self.accepted, len(self.results)) + \
               MoveChange._indent('\n'.join(map(str, self.subchanges)))


class NonCanonicalConditionalSequentialMoveChange(
                                            ConditionalSequentialMoveChange):
    """ Special move change for reactive flux simulation.

    This move change inherits from :class:`.ConditionalSequentialMoveChange`
    and returns the outcome of the last subchange.
    """
    @property
    def canonical(self):
        return self.subchanges[-1]


class SubMoveChange(MoveChange):
    """
    A helper MoveChange that represents the application of a submover.

    The raw implementation delegates all to the subchange
    """
    def __init__(self, subchange, mover=None, details=None):
        """
        Parameters
        ----------
        subchange : MoveChange
            the actual subchange used by this wrapper PMC
        mover
        details

        Attributes
        ----------
        subchange
        mover
        details
        """
        super(SubMoveChange, self).__init__(mover=mover, details=details)
        self.subchanges = [subchange]

    def _get_results(self):
        return self.subchange.results

    def _get_trials(self):
        return self.subchange.trials

    def __str__(self):
        # Defaults to use the name of the used mover
        return self.mover.__class__.__name__[:-5] + ' :\n' + MoveChange._indent(str(self.subchange))


class RandomChoiceMoveChange(SubMoveChange):
    """
    A MoveChange that represents the application of a mover chosen randomly
    """

    # This class is empty since all of the decision is specified by the mover
    # and it requires no additional logic to decide if it is accepted.

class FilterByEnsembleMoveChange(SubMoveChange):
    """
    A MoveChange that filters out all samples not in specified ensembles
    """

    # TODO: Question: filter out also trials not in the ensembles? I think so,
    # because we are only interested in trials that could be relevant, right?

    def _get_results(self):
        all_samples = self.subchange.results

        filtered_samples = list(filter(
            lambda s: s.ensemble in self.mover.ensembles,
            all_samples
        ))

        return filtered_samples

    def _get_trials(self):
        all_samples = self.subchange.trials

        filtered_samples = list(filter(
            lambda s: s.ensemble in self.mover.ensembles,
            all_samples
        ))

        return filtered_samples


    def __str__(self):
        return 'FilterMove : allow only ensembles [%s] from sub moves : %s : %d samples\n' % \
               (str(self.mover.ensembles), self.accepted, len(self.results)) + \
               MoveChange._indent(str(self.subchange))



class FilterSamplesMoveChange(SubMoveChange):
    """
    A MoveChange that keeps a selection of the underlying samples
    """

    def _get_results(self):
        sample_set = self.subchange.results

        # allow for negative indices to be picked, e.g. -1 is the last sample
        samples = [idx % len(sample_set) for idx in self.mover.selected_samples]

        return samples

    def __str__(self):
        return 'FilterMove : pick samples [%s] from sub moves : %s : %d samples\n' % \
               (str(self.mover.selected_samples), self.accepted, len(self.results)) + \
               MoveChange._indent(str(self.subchange))


class KeepLastSampleMoveChange(SubMoveChange):
    """
    A MoveChange that only keeps the last generated sample.

    This is different from using `.reduced` which will only change the
    level of detail that is stored. This MoveChange will actually remove
    potential relevant samples and thus affect the outcome of the new
    SampleSet. To really remove samples also from storage you can use
    this MoveChange in combination with `.closed` or `.reduced`

    Notes
    -----
    Does the same as `FilterSamplesMoveChange(subchange, [-1], False)`

    I think we should try to not use this. It would be better to make submoves
    and finally filter by relevant ensembles. Much like running a function
    with local variables/local ensembles.
    """

    def _get_results(self):
        samples = self.subchange.results
        if len(samples) > 1:
            samples = [samples[-1]]

        return samples

    def __str__(self):
        return 'Restrict to last sample : %s : %d samples\n' % \
               (self.accepted, len(self.results)) + \
               MoveChange._indent(str(self.subchange))


class PathSimulatorMoveChange(SubMoveChange):
    """
    A MoveChange that just wraps a subchange and references a PathSimulator
    """

    def __str__(self):
        return 'PathSimulatorStep : %s : Step # %d with %d samples\n' % \
               (str(self.mover.pathsimulator.cls), self.details.step, len(self.results)) + \
               MoveChange._indent(str(self.subchange))
