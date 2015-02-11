__author__ = 'jan-hendrikprinz'

import openpathsampling as paths
from openpathsampling.todict import restores_as_full_object

@restores_as_full_object
class MovePath(object):
    '''
    MovePath is essentially a list of samples, with a few conveniences.  It
    can be treated as a list of samples (using, e.g., .append), or as a
    dictionary of ensembles mapping to a list of samples, or as a dictionary
    of replica IDs to samples. Any type is allowed as a replica ID except
    Sample or Ensemble.

    The dictionaries ensemble_dict and replica_dict are conveniences which
    should be kept consistent by any method which modifies the container.
    They do not need to be stored.

    Notes
    -----
    Current implementation is as an unordered set. Therefore we don't
    have some of the convenient tools in Python sequences (e.g.,
    slices). On the other hand, I'm not sure whether that is meaningful
    here.

    Attributes
    ----------
    samples : list of Sample
        The samples included in this set.
    ensemble_dict : dict
        A dictionary with Ensemble objects as keys and lists of Samples as
        values.
    replica_dict : dict
        A dictionary with replica IDs as keys and lists of Samples as values
    '''

    @staticmethod
    def _indent(s):
        spl = s.split('\n')
        spl = [' |  ' + p if p[0] == ' ' else ' +- ' + p for p in spl]
        return '\n'.join(spl)

    def __init__(self, accepted=True, mover=None):
        self._accepted = accepted
        self._destination = None
        self._local_samples = []
        self._collapsed_samples = None
        self._samples = None
        self.mover = mover

    def to_dict(self):
        return {
            'accepted' : self.accepted,
            'mover' : self.mover,
        }

    @property
    def opened(self):
        """
        Return the full MovePath object

        Notes
        -----
        The full Movepath can only be returned if it has been generated and
        is still in memory. A collapsed Movepath that is loaded does NOT
        contain the full information anymore. This is the whole purpose of
        only storing essential information.
        """
        return self

    @property
    def closed(self):
        """
        Return a closed MovePath object copy

        This will return a MovePath object that will still know the used
        mover and all relevant samples. All underlying information about the
        move is hidden and will not be stored
        """
        obj = CollapsedMovePath(samples=self.collapsed_samples, mover=self.mover)
        obj._movepath = self

        return obj

    def reduced(self, selected_samples = None, use_all_samples=False):
        """
        Reduce the underlying MovePath to a subset of relevant samples

        This can be used to reduce a MovePath (and everything below) to
        a movepath that only uses specific samples that can be picked
        from the list of all returned samples. Note that this can be
        problematic since a move might not return a different number of
        samples each time it is run.

        Parameters
        ----------
        selected_samples : list of int
            list of integer indices to be kept from the list of returned
            samples in this movepath. Slicing is not allowed, but negative
            indices are and conforms to the usual python convention (-1
            is the last samples, etc.)
        use_all_samples : bool
            if `True` the selected samples will be chosen from the list of
            all created samples (accepted and rejected ones), otherwise only
            the accepted ones will be chosen from. This is the default and
            corresponds to chose from `.samples`
        """

        # @DWHS do we want to allow also to select rejected samples? Meaning,
        # that we could allow the user to pick samples from .all_samples
        # or .samples

        # the check for collapsed_samples makes sure that at least the
        # relevant ones are present

        if selected_samples is None:
            # this case corresponds to .closed
            return self.closed

        if use_all_samples:
            # choose all generated samples
            sample_set = self.all_samples
        else:
            # chose only accepted ones!
            sample_set = self.samples

        # allow for negative indices to be picked, e.g. -1 is the last sample
        if len(selected_samples) > 0:
            selected_samples = [ idx % len(sample_set) for idx in selected_samples]

        samples = [
            samp for idx, samp in enumerate(sample_set)
            if idx in selected_samples or samp in self.collapsed_samples
        ]

        obj = CollapsedMovePath(samples = samples, mover=self.mover)
        obj._movepath = self

        return obj

    @property
    def collapsed_samples(self):
        """
        Return a collapsed set of samples with non used samples removed

        This is the minimum required set of samples to keep the `MovePath`
        correct and allow to target sample set to be correctly created.
        These are the samples used by `.closed`

        Example
        -------
        Assume that you run 3 shooting moves for replica #1. Then only the
        last of the three really matters for the target sample_set since #1
        will be replaced by #2 which will be replaced by #3. So this function
        will return only the last sample.

        See also
        --------
        MovePath.closed
        MovePath.reduced()

        """
        if self._collapsed_samples is None:
            s = paths.SampleSet([]).apply_samples(self.samples)

            # keep order just for being thorough
            self._collapsed_samples = [
                samp for samp in self.samples
                if samp in s
            ]

        return self._collapsed_samples

    @property
    def local_samples(self):
        """
        A list of the samples that are needed to update the new sampleset

        Notes
        -----
        These samples are purely for this MovePath node and not for any
        underlying nodes. It is effectively only used by SampleMovePath
        """
        return self._local_samples

    @property
    def all_samples(self):
        """
        Returns a list of all samples generated during the PathMove.

        This includes all accepted and rejected samples (which does NOT
        include hidden samples yet)

        TODO: Decide, if we want hidden samples to be included or not
        """
        return self._local_samples

    @property
    def accepted(self):
        """
        Returns if this particular move was accepted.

        Mainly used for rejected samples.
        """
        if self._accepted is None:
            self._accepted = self._get_accepted()

        return self._accepted

    def _get_accepted(self):
        """
        The function that determines if a MovePath was accepted.

        Is overridden by SequentialConditionalMovePath, etc.
        """
        return True

    @property
    def __add__(self, other):
        return SequentialMovePath([self, other])

    def apply_to(self, other):
        """
        Standard apply is to apply the list of samples contained
        """
        new_sampleset = paths.SampleSet(other).apply_samples(self.samples)
        new_sampleset.movepath = self

    @property
    def samples(self):
        """
        Returns a list of all samples that are accepted in this move

        This contains unnecessary, but accepted samples, too.
        """
        if self._samples is None:
            self._samples = self._get_samples()

        return self._samples

    def _get_samples(self):
        """
        Determines all relevant accepted samples for this move

        Includes all accepted samples also from submoves
        """
        if self.accepted:
            return self._local_samples
        else:
            return []

    def __iter__(self):
        """
        Allow to iterate over all accepted samples in the move

        This effectively iterates over `self.samples`
        """
        for sample in self.samples:
            yield sample

    def __len__(self):
        """
        Returns the number of contained accepted samples

        Shortcut for `len(self.samples)`
        """
        return len(self.samples)

    def __contains__(self, item):
        """
        Check, if a particular sample is among the accepted samples
        """
        return (item in self.samples)

    def __str__(self):
        if self.accepted:
            return 'SampleMove : %s : %s : %d samples' % (self.mover.cls, self.accepted, len(self._local_samples)) + ' ' + str(self._local_samples) + ''
        else:
            return 'SampleMove : %s : %s :[]' % (self.mover.cls, self.accepted)


@restores_as_full_object
class EmptyMovePath(MovePath):
    """
    A MovePath representing no changes
    """
    def __init__(self, mover=None):
        super(EmptyMovePath, self).__init__(accepted=True, mover=mover)

    def apply_to(self, other):
        return other

    def __str__(self):
        return ''

    def to_dict(self):
        return {}

    @property
    def all_samples(self):
        return []



@restores_as_full_object
class SampleMovePath(MovePath):
    """
    A MovePath representing the application of samples.

    This is the most common MovePath and all other moves use this
    as leaves and on the lowest level consist only of `SampleMovePath`
    """
    def __init__(self, samples, accepted=True, mover=None):
        super(SampleMovePath, self).__init__(accepted=accepted, mover=mover)
        if samples is None:
            return

        if type(samples) is paths.Sample:
            samples = [samples]


        self._local_samples.extend(samples)

    def to_dict(self):
        return {
            'accepted' : self.accepted,
            'mover' : self.mover,
            'samples' : self._local_samples
        }

    def apply_to(self, other):
        """
        Standard apply is to apply the list of samples contained
        """
        return paths.SampleSet(other).apply_samples(self._local_samples)


@restores_as_full_object
class CollapsedMovePath(SampleMovePath):
    """
    Represent a collapsed MovePath that has potential hidden sub moves
    """
    @property
    def opened(self):
        if hasattr(self, '_movepath') and self._movepath is not None:
            return self._movepath
        else:
            return self

    def closed(self):
        return self

    @property
    def collapsed_samples(self):
        if self._collapsed_samples is None:

            self._collapsed_samples = self._local_samples

        return self._collapsed_samples

    def __str__(self):
        if self.mover is not None:
            return '%s [collapsed] : %d samples' % (self.mover.cls, len(self._local_samples)) + ' ' + str(self._local_samples) + ''
        else:
            return '%s [collapsed] : %d samples' % ('CollapsedMove', len(self._local_samples)) + ' ' + str(self._local_samples) + ''


@restores_as_full_object
class RandomChoiceMovePath(MovePath):
    """
    A MovePath that represents the application of a mover chosen randomly
    """
    def __init__(self, movepath, mover=None):
        super(RandomChoiceMovePath, self).__init__(mover=mover)
        self.movepath = movepath

    def to_dict(self):
        return {
            'movepath' : self.movepath
        }

    def _get_samples(self):
        return self.movepath.samples

    def apply_to(self, other):
        return self.movepath.apply_to(other)

    def __str__(self):
        return 'RandomChoice :\n' + MovePath._indent(str(self.movepath))

    @property
    def all_samples(self):
        return self.movepath.all_samples


@restores_as_full_object
class SequentialMovePath(MovePath):
    """
    SequentialMovePath has no own samples, only inferred Sampled from the
    underlying MovePaths
    """
    def __init__(self, movepaths, mover=None):
        super(SequentialMovePath, self).__init__(accepted=None, mover=mover)
        self.movepaths = movepaths

    def to_dict(self):
        return {
            'movepath' : self.movepaths
        }

    def _get_samples(self):
        samples = []
        for movepath in self.movepaths:
            samples = samples + movepath.samples
        return samples

    @property
    def all_samples(self):
        samples = []
        for movepath in self.movepaths:
            samples = samples + movepath.all_samples
        return samples


    def apply_to(self, other):
        sampleset = other

        for movepath in self.movepaths:
            sampleset = movepath.apply_to(sampleset)

        return sampleset

    def __str__(self):
        return 'SequentialMove : %s : %d samples\n' % \
               (self.accepted, len(self.samples)) + \
               MovePath._indent('\n'.join(map(str, self.movepaths)))

@restores_as_full_object
class PartialAcceptanceSequentialMovePath(SequentialMovePath):
    """
    PartialAcceptanceSequentialMovePath has no own samples, only inferred
    Sampled from the underlying MovePaths
    """

    def _get_samples(self):
        changes = []
        for movepath in self.movepaths:
            if movepath.accepted:
                changes.extend(movepath.samples)
            else:
                break

        return changes

    def apply_to(self, other):
        sampleset = other

        for movepath in self.movepaths:
            if movepath.accepted:
                sampleset = movepath.apply_to(sampleset)
            else:
                break

        return sampleset

    def __str__(self):
        return 'PartialAcceptanceMove : %s : %d samples\n' % \
               (self.accepted, len(self.samples)) + \
               MovePath._indent('\n'.join(map(str, self.movepaths)))

@restores_as_full_object
class ConditionalSequentialMovePath(SequentialMovePath):
    """
    ConditionalSequentialMovePath has no own samples, only inferred Samples
    from the underlying MovePaths
    """

    def apply_to(self, other):
        sampleset = other

        for movepath in self.movepaths:
            if movepath.accepted:
                sampleset = movepath.apply_to(sampleset)
            else:
                return other

        return sampleset

    def _get_samples(self):
        changes = []
        for movepath in self.movepaths:
            if movepath.accepted:
                changes.extend(movepath)
            else:
                return []

        return changes

    def _get_accepted(self):
        for movepath in self.movepaths:
            if not movepath.accepted:
                return False

        return True

    def __str__(self):
        return 'ConditionalSequentialMove : %s : %d samples\n' % \
               (self.accepted, len(self.samples)) + \
               MovePath._indent( '\n'.join(map(str, self.movepaths)))


@restores_as_full_object
class FilterSamplesMovePath(MovePath):
    """
    A MovePath that keeps a selection of the underlying samples
    """

    def __init__(self, movepath, selected_samples, use_all_samples=False, mover=None):
        super(KeepLastSampleMovePath, self).__init__(mover=mover)
        self.movepath = movepath
        self.selected_samples = selected_samples
        self.use_all_samples = use_all_samples


    def _get_samples(self):
        if self.use_all_samples:
            # choose all generated samples
            sample_set = self.all_samples
        else:
            # chose only accepted ones!
            sample_set = self.samples

        # allow for negative indices to be picked, e.g. -1 is the last sample
        samples = [ idx % len(sample_set) for idx in self.selected_samples]

        return samples

    def __str__(self):
        return 'FilterMove : pick samples [%s] from sub moves : %s : %d samples\n' % \
               (str(self.selected_samples), self.accepted, len(self.samples)) + \
               MovePath._indent( str(self.movepath) )

    def to_dict(self):
        return {
            'movepath' : self.movepath,
            'selected_samples' : self.selected_samples,
            'use_all_samples' : self.use_all_samples
        }

@restores_as_full_object
class KeepLastSampleMovePath(MovePath):
    """
    A MovePath that only keeps the last generated sample.

    This is different from using `.reduced` which will only change the
    level of detail that is stored. This MovePath will actually remove
    potential relevant samples and thus affect the outcome of the new
    SampleSet. To really remove samples also from storage you can use
    this MovePath in combination with `.closed` or `.reduced`

    Notes
    -----
    Does the same as `FilterSamplesMovePath(movepath, [-1], False)`
    """
    def __init__(self, movepath, mover=None):
        super(KeepLastSampleMovePath, self).__init__(mover=mover)
        self.movepath = movepath

    def _get_samples(self):
        samples = self.movepath.samples
        if len(samples) > 1:
            samples = [samples[-1]]

        return samples

    def __str__(self):
        return 'Restrict to last sample : %s : %d samples\n' % \
               (self.accepted, len(self.samples)) + \
               MovePath._indent( str(self.movepath) )

    def to_dict(self):
        return {
            'movepath' : self.movepath
        }

@restores_as_full_object
class CalculationMovePath(MovePath):
    """
    A MovePath that just wraps a movepath and references a Calculation
    """

    def __init__(self, movepath, calculation=None, step=-1):
        super(CalculationMovePath, self).__init__(mover=movepath.mover)
        self.movepath = movepath
        self.calculation = calculation
        self.step = step

    def _get_samples(self):
        return self.movepath.samples

    def __str__(self):
        return 'CalculationStep : %s : Step # %d with %d samples\n' % \
               (str(self.calulation.cls), self.step, len(self.samples)) + \
               MovePath._indent( str(self.movepath) )

    def to_dict(self):
        return {
            'movepath' : self.movepath,
            'calculation' : self.calculation,
            'step' : self.step
        }