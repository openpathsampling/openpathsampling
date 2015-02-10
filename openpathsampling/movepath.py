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
        self._samples = None
        self.mover = mover

    def to_dict(self):
        return {
            'accepted' : self.accepted,
            'mover' : self.mover,
        }

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
            return 'SampleMove : %s : %s : %d samples' % (self.mover.__class__.__name__, self.accepted, len(self._local_samples)) + ' ' + str(self._local_samples) + ''
        else:
            return 'SampleMove : %s : %s :[]' % (self.mover.__class__.__name__, self.accepted)


@restores_as_full_object
class EmptyMovePath(MovePath):
    def __init__(self):
        super(EmptyMovePath, self).__init__(accepted=True)

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
class RandomChoiceMovePath(MovePath):
    """
    RandomMoveMovePath contains only a reference to the underlying used
    MovePath
    """
    def __init__(self, movepath):
        super(RandomChoiceMovePath, self).__init__()
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
    def __init__(self, movepaths):
        super(SequentialMovePath, self).__init__(accepted=None)
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
class PartialMovePath(SequentialMovePath):
    """
    PartialMovePath has no own samples, only inferred Sampled from the
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
class ExclusiveMovePath(SequentialMovePath):
    """
    ExclusiveMovePath has no own samples, only inferred Sampled from the
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

@restores_as_full_object
class KeepLastSampleMovePath(MovePath):
    def __init__(self, movepath):
        super(KeepLastSampleMovePath, self).__init__()
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


class MovePathIdea(object):
    """
    Information Class to contain the actual route/path of moves taken in
    a complex (Path)Mover object

    The contains all information necessary to
    (a) turn the initial_sampleset into the final_sampleset
    (b) contains the result of all submovers used in the generation of
        the final MovePath, this includes the order in which the submoves
        are called, if their result is the same as before (thus accepted)
        and what samples were apply to the current (intermediate) MovePath

    - The intermediate samplesets are not stored but used internally to track
    results.
    - Each submove has a MovePath as result and samplesets are changed using
    samples
    -

    Examples
    --------
    >>> mover = Sequential([moveA1,
    >>>     Progressive([moveB1, moveB2]),
    >>>     Exclusive([moveC1, moveC2]),
    >>>     moveA2
    >>> ])
    >>> bucket2 = mover(bucket1)
    >>> movepath = bucket2.movepath
    >>> print movepath

    # All attempted move possibilities

    1: [moveA1:y, moveB1:y, moveB2:y, moveC1:y, moveC2:y, moveA2:y]
    2: [moveA1:y, moveB1:y, moveB2:y, moveC1:y, moveC2:n, moveA2:y]
    3: [moveA1:y, moveB1:y, moveB2:y, moveC1:n,           moveA2:y]
    4: [moveA1:y, moveB1:y, moveB2:n, moveC1:y, moveC2:y, moveA2:y]
    5: [moveA1:y, moveB1:y, moveB2:n, moveC1:y, moveC2:n, moveA2:y]
    6: [moveA1:y, moveB1:y, moveB2:n, moveC1:n,           moveA2:y]
    7: [moveA1:y, moveB1:n,           moveC1:y, moveC2:y, moveA2:y]
    8: [moveA1:y, moveB1:n,           moveC1:y, moveC2:n, moveA2:y]
    9: [moveA1:y, moveB1:n,           moveC1:n,           moveA2:y]

    # Accepted Move possibilities
    1: [moveA1:y, moveB1:y, moveB2:y, moveC1:y, moveC2:y, moveA2:y]
    2: [moveA1:y, moveB1:y, moveB2:y, moveC1:y,           moveA2:y]
    3: [moveA1:y, moveB1:y, moveB2:y,                     moveA2:y]
    4: [moveA1:y, moveB1:y,           moveC1:y, moveC2:y, moveA2:y]
    5: [moveA1:y, moveB1:y,           moveC1:y,           moveA2:y]
    6: [moveA1:y, moveB1:y,                               moveA2:y]
    7: [moveA1:y,                     moveC1:y, moveC2:y, moveA2:y]
    8: [moveA1:y,                     moveC1:y,           moveA2:y]
    9: [moveA1:y,                                         moveA2:y]

    # All possible move paths
    1: S([moveA1:y, P([moveB1:y, moveB2:y]):y, E([moveC1:y, moveC2:y]):y, moveA2:y]):y
    2: S([moveA1:y, P([moveB1:y, moveB2:y]):y, E([moveC1:y, moveC2:n]):n, moveA2:y]):y
    3: S([moveA1:y, P([moveB1:y, moveB2:y]):y, E([moveC1:n]):n,           moveA2:y]):y
    4: S([moveA1:y, P([moveB1:y, moveB2:n]):y, E([moveC1:y, moveC2:y]):y, moveA2:y]):y
    5: S([moveA1:y, P([moveB1:y, moveB2:n]):y, E([moveC1:y, moveC2:n]):n, moveA2:y]):y
    6: S([moveA1:y, P([moveB1:y, moveB2:n]):y, E([moveC1:n]):n,           moveA2:y]):y
    7: S([moveA1:y, P([moveB1:n]):y,           E([moveC1:y, moveC2:y]):y, moveA2:y]):y
    8: S([moveA1:y, P([moveB1:n]):y,           E([moveC1:y, moveC2:n]):n, moveA2:y]):y
    9: S([moveA1:y, P([moveB1:n]):y,           E([moveC1:n]):n,           moveA2:y]):y

    # Number of all moves
    1: 6
    2: 6
    3: 5
    4: 6
    5: 6
    6: 5
    7: 5
    8: 5
    9: 4

    # Number of accepted moves
    1: 6
    2: 5
    3: 4
    4: 5
    5: 4
    6: 3
    7: 4
    8: 3
    9: 2

    # Examples for no #8

    S([
        moveA1:y:[sample1, sample2],
        P([
            moveB1:y:[sample3, sample4],
            moveB2:n:[sample5, sample6]
        ]):y:[sample3, sample4],
        E([
            moveC1:y:[sample7, sample8, sample9],
            moveC2:n:[sample10, sample11]
        ]):n:[],
        moveA2:y:[sample12]
    ]):y:[sample1, sample2, sample3, sample4, sample12]

    MovePaths form a group like structure with our known concatenation of
    applying all samples in a row. But this requires that application of
    Samples is unique. Can we assure that? This requires to set an initial
    sample that is to be replaced by another or nothing. Adding just means to
    start from an emtpy sample. The randomness need to be stored in the sample.
    The we can really do

    setA * setB = setC

    setA * setB * setC = setA * (setB * setC) = (setA * setB) * setC

    setA * setB = []

    We need MovePath to know the original MovePath and a MovePath
    that contain the moves, which can be nested. So A Move takes the full
    MoveSample and Returns the object with added Moves which are itself
    a MovePath.

    SeqentialMovePath(
        initial_set, [
        MovePath(
            moveA1:y:initial_set + [sample1, sample2]
        ),
        PartialMovePath(
        [
            moveB1:y:initial_set + [sample1, sample2] + [sample3, sample4],
            moveB2:n:initial_set + [sample1, sample2] + [sample3, sample4] + [sample5, sample6]
        ]):y:initial_set + [sample1, sample2] + [sample3, sample4],
        ExcluiveMovePath([
            moveC1:y:initial_set + [sample1, sample2] + [sample3, sample4] + [sample7, sample8, sample9],
            moveC2:n:initial_set + [sample1, sample2] + [sample3, sample4] + [sample7, sample8, sample9] + [sample10, sample11]
        ]):n:[],
        moveA2:y:initial_set + [sample1, sample2] + [sample3, sample4] + [sample12]
    ]):y:initial_set + [sample1, sample2] + [sample3, sample4] + [sample12]


    Attributes
    ----------
    initial_sampleset : MovePath
        the initial MovePath where the path originates
    final_sampleset : MovePath
        the final MovePath where the path ends
    pathmover : PathMover
        a reference to the full PathMover object called that generated
        this path

    """

    def __init__(self):
        self.initial_sampleset = None
        self.final_sampleset = None
        self.pathmover = None
        self.moves = None
        pass

    def samples(self):
        """
        Returns all samples generated along the movepath

        :return:
        """

        list_of_samples = []

        return list_of_samples

    def relevent_samples(self):
        """
        Returns all samples that are relevant in the generation of the
        new MovePath from the old one. Overwritten samples are removed.

        Usually this set should be the same as samples

        :return:
        """

        list_of_samples = []

        return list_of_samples

    def count_movers(self, mover):
        """
        Return the number of times a specific mover has been called.

        :return:
        """

        counts = 0

        return counts


    def count_accepts(self, mover):
        """
        Return number of times a specific mover has been accepted.

        :param mover:
        :return:
        """

        counts = 0

        return counts

    def count_rejects(self, mover):
        """
        Return number of times a specific mover has been rejected.

        :param mover:
        :return:
        """

        counts = 0

        return counts
