__author__ = 'jan-hendrikprinz'

import openpathsampling as paths

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

    def __init__(self, origin = None, accepted=True):
        self._accepted = accepted
        self._destination = None
        self._origin = None
        self._samples = []

    @property
    def samples(self):
        """
        A list of the contained samples for this particular move
        """
        return self._samples

    @property
    def accepted(self):
        if self._accepted is None:
            self._accepted = self._get_accepted()

        return self._accepted

    def _get_accepted(self):
        return True

    @property
    def __radd__(self, other):
        if other is self.origin:
            if self._destination is None:
                self._destination = self._apply(other)
            return self._destination
        else:
            return other

    def _apply(self, other):
        """
        Standard apply is to apply the list of samples contained
        """
        return paths.SampleSet(other).apply_samples(self._samples)

    @property
    def changes(self):
        """
        Returns a list of all samples that should be applied to
        """
        if self._changes is None:
            self._get_changes()

        return self._changes

    def _get_changes(self):
        """
        A normal movepath
        """
        self._changes = self._samples

    def __iter__(self):
        for sample in self.changes:
            yield sample

    def __len__(self):
        return len(self.changes)

    def __contains__(self, item):
        """
        Check, if a particular
        """
        return (item in self.changes)


class EmptyMovePath(MovePath):
    def __init__(self):
        super(EmptyMovePath, self).__init__(accepted=True)

    def _apply(self, other):
        return other


class SampleMovePath(MovePath):
    """
    Sample Move Path

    Most common
    """
    def __init__(self, origin, samples, accepted=True):
        super(SampleMovePath, self).__init__(origin=origin, accepted=accepted)
        if samples is None:
            return

        if type(samples) is paths.Sample:
            samples = [samples]

        self._samples.extend(samples)

    def _apply(self, other):
        """
        Standard apply is to apply the list of samples contained
        """
        return paths.SampleSet(other).apply_samples(self._samples)



class RandomChoiceMovePath(MovePath):
    """
    RandomMoveMovePath contains only a reference to the underlying used
    MovePath
    """
    def __init__(self, origin, movepath):
        super(RandomChoiceMovePath, self).__init__(origin=origin)
        self.movepath = movepath

    def _get_changes(self):
        return self.movepath._get_changes()

    def _apply(self, other):
        return self.movepath._apply(other)


class SequentialMovePath(MovePath):
    """
    SequentialMovePath has no own samples, only inferred Sampled from the
    underlying MovePaths
    """
    def __init__(self, origin, movepaths):
        super(SequentialMovePath, self).__init__(origin=origin, accepted=None)
        self.movepaths = movepaths

    def _get_changes(self):
        changes = []
        map(changes.extend, self.movepaths)
        return changes

    def _apply(self, other):
        sampleset = other

        for movepath in self.movepaths:
            sampleset = movepath._apply(sampleset)

        return sampleset


class PartialMovePath(SequentialMovePath):
    """
    PartialMovePath has no own samples, only inferred Sampled from the
    underlying MovePaths
    """

    def _get_changes(self):
        changes = []
        for sampleset in self.movepaths:
            if sampleset.accepted:
                changes.extend(sampleset)
            else:
                break

        return changes

    def _apply(self, other):
        sampleset = other

        for movepath in self.movepaths:
            if movepath.accepted:
                sampleset = movepath._apply(sampleset)
            else:
                break

        return sampleset

class ExclusiveMovePath(SequentialMovePath):
    """
    ExclusiveMovePath has no own samples, only inferred Sampled from the
    underlying MovePaths
    """

    def _apply(self, other):
        sampleset = other

        for movepath in self.movepaths:
            if movepath.accepted:
                sampleset = movepath._apply(sampleset)
            else:
                return other

        return sampleset

    def _get_changes(self):
        changes = []
        for sampleset in self.samplesets:
            if sampleset.accepted:
                changes.extend(sampleset)
            else:
                return []

        return changes

    def _get_accepted(self):
        for sampleset in self.samplesets:
            if not sampleset.accepted:
                return False

        return True



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
