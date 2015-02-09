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

    Note
    ----
        Current implementation is as an unordered set. Therefore we don't
        have some of the convenient tools in Python sequences (e.g.,
        slices). On the other hand, I'm not sure whether that is meaningful
        here.

        JHP: I changed it to be immutable! This makes it slower but better
        tractable. Just a try

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
        """
        return self

    @property
    def closed(self):
        """
        Return a closed MovePath object copy
        """
        obj = CollapsedMovePath(samples=self.collapsed_samples, mover=self.mover)
        obj._movepath = self
        return obj

    @property
    def collapsed_samples(self):
        """
        Return a collapsed set of samples with non used samples removed
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
        """
        return self._local_samples

    @property
    def all_samples(self):
        """
        Returns a list of all samples generated during the PathMove.

        This includes all rejected samples
        """
        return self._local_samples

    @property
    def accepted(self):
        if self._accepted is None:
            self._accepted = self._get_accepted()

        return self._accepted

    def _get_accepted(self):
        return True

    @property
    def __add__(self, other):
        return SequentialMovePath([self, other])

    def apply_to(self, other):
        """
        Standard apply is to apply the list of samples contained
        """
        return paths.SampleSet(other).apply_samples(self._local_samples)

    @property
    def samples(self):
        """
        Returns a list of all samples that should be applied to
        """
        if self._samples is None:
            self._samples = self._get_samples()

        return self._samples

    def _get_samples(self):
        """
        A normal movepath
        """
        if self.accepted:
            return self._local_samples
        else:
            return []

    def __iter__(self):
        for sample in self.samples:
            yield sample

    def __len__(self):
        return len(self.samples)

    def __contains__(self, item):
        """
        Check, if a particular
        """
        return (item in self.samples)

    def __str__(self):
        if self.accepted:
            return 'SampleMove : %s : %s : %d samples' % (self.mover.cls, self.accepted, len(self._local_samples)) + ' ' + str(self._local_samples) + ''
        else:
            return 'SampleMove : %s : %s :[]' % (self.mover.cls, self.accepted)


@restores_as_full_object
class EmptyMovePath(MovePath):
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
    Sample Move Path

    Most common
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
        """
        Return a collapsed set of samples with non used samples removed
        """
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
    RandomMoveMovePath contains only a reference to the underlying used
    MovePath
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

        if self.accepted:
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
            'movepaths' : self.movepaths
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
    PartialAcceptanceSequentialMovePath has no own samples, only inferred Sampled from the
    underlying MovePaths
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
        return 'PartialMove : %s : %d samples\n' % \
               (self.accepted, len(self.samples)) + \
               MovePath._indent('\n'.join(map(str, self.movepaths)))

@restores_as_full_object
class ConditionalSequentialMovePath(SequentialMovePath):
    """
    ConditionalSequentialMovePath has no own samples, only inferred Sampled from the
    underlying MovePaths
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
        return 'ExclusiveMove : %s : %d samples\n' % \
               (self.accepted, len(self.samples)) + \
               MovePath._indent( '\n'.join(map(str, self.movepaths)))

class KeepLastSampleMovePath(MovePath):
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
