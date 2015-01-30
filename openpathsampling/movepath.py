__author__ = 'jan-hendrikprinz'

class SampleSet(object):
    '''
    SampleSet is essentially a list of samples, with a few conveniences.  It
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

    # TODO: Can a sample be several times in an SampleSet? Why not?

    def __init__(self, samples = None, predecessor = None, accepted=True):
        self.predecessor = predecessor
        self._accepted = accepted
        self._ensemble_dict = None
        self._replica_dict = None
        self._samples = []
        self._destination = None
        self._origin = None

        if samples is None:
            return

        if type(samples) is Sample:
            samples = [samples]

        self._samples.extend(samples)

    @property
    def samples(self):
        if self.predecessor is None:
            return self._samples
        else:
            return self.predecessor + self._samples

    @property
    def accepted(self):
        if self._accepted is None:
            self._accepted = self._get_accepted()

        return self._accepted

    def _get_accepted(self):
        return True

    @property
    def destination(self):
        if self._destination is None:
            self._destination = self._get_destination()

        return self._destination

    def _get_destination(self):
        return self.predecessor + SampleSet(self._samples, self)

    @property
    def origin(self):
        if self._origin is None:
            self._origin = self._get_origin()

        return self._origin


    def _get_origin(self):
        return self.predecessor.destination

    @property
    def changes(self):
        if self._changes is None:
            self._get_changes()

        return self._destination

    def _get_changes(self):
        self._changes = self._samples


    @property
    def ensemble_dict(self):
        if self._ensemble_dict is None:
            self._ensemble_dict = self._get_ensemble_dict()

        return self._ensemble_dict

    def _get_ensemble_dict(self):
        """
        Returns the dictionary of ensembles and their samples but not cached
        :return:
        """
        ensembles = set([sample.ensemble for sample in self.samples])
        print ensembles
        return { sample.ensemble : [sample for sample in self.samples if sample.ensemble is ensemble] for ensemble in ensembles}


    @property
    def replica_dict(self):
        if self._replica_dict is None:
            self._replica_dict = self._get_replica_dict()

        return self._replica_dict

    def _get_replica_dict(self):
        """
        Returns the dictionary of replica and their samples but not cached
        :return:
        """
        replicas = set([sample.replica for sample in self.samples])
        return { sample.replica : [sample for sample in self.samples if sample.replica is replica] for replica in replicas}

    def __plus__(self, other):
        if other.predecessor is self:
            newset = self.copy()
            for sample in other._samples:
                if sample not in self._samples:
                    self._append(sample)

            return newset
        else:
            raise ValueError('Incompatible SampleSets')

    def __getitem__(self, key):
        if isinstance(key, paths.Ensemble):
            print self.ensemble_dict
            return random.choice(self.ensemble_dict[key])
        else:
            print self.replica_dict
            return random.choice(self.replica_dict[key])

    def __iter__(self):
        for sample in self.samples:
            yield sample

    def __len__(self):
        return len(self.samples)

    def __contains__(self, item):
        return (item in self.samples)

    def all_from_ensemble(self, ensemble):
        try:
            return self.ensemble_dict[ensemble]
        except KeyError:
            return []

    def all_from_replica(self, replica):
        try:
            return self.replica_dict[replica]
        except KeyError:
            return []

    def replica_list(self):
        '''Returns the list of replicas IDs in this SampleSet'''
        return self.replica_dict.keys()

    def ensemble_list(self):
        '''Returns the list of ensembles in this SampleSet'''
        return self.ensemble_dict.keys()


class EmptySampleSet(SampleSet):
    def __init__(self):
        super(EmptySampleSet, self).__init__(accepted=True)

    def _get_changes(self):
        return []

    def _get_origin(self):
        return []

    def _get_final(self):
        return []

class InitialSampleSet(SampleSet):
    """
    A SampleSet with empty initial referenced sampleset
    """
    def __init__(self, samples=None):
        super(InitialSampleSet, self).__init__(samples=samples, accepted=True)

    def _get_changes(self):
        return []

    def _get_origin(self):
        return []

    def _get_final(self):
        return []

class RandomMoveSampleSet(SampleSet):
    """
    RandomMoveSampleSet contains only a reference to the underlying used
    SampleSet
    """
    def __init__(self, sampleset, predecessor = None):
        super(RandomMoveSampleSet, self).__init__(predecessor=predecessor)
        self.sampleset = sampleset

    def _get_changes(self):
        return self.sampleset._get_changes

    def _get_origin(self):
        return self.sampleset._get_origin

    def _get_changes(self):
        return self.sampleset._get_changes

    def _get_changes(self):
        return self.sampleset._get_changes


class SequentialSampleSet(SampleSet):
    """
    SequentialSampleSet has no own samples, only inferred Sampled from the
    underlying SampleSets
    """
    def __init__(self, samplesets, predecessor = None):
        super(SequentialSampleSet, self).__init__(predecessor=predecessor)
        self.samplesets = samplesets

    def _get_changes(self):
        changes = []
        map(changes.extend, self.samplesets)
        return changes

class PartialSampleSet(SequentialSampleSet):
    """
    PartialSampleSet has no own samples, only inferred Sampled from the
    underlying SampleSets
    """

    def _get_changes(self):
        changes = []
        for sampleset in self.samplesets:
            if sampleset.accepted:
                changes.extend(sampleset)
            else:
                break

        return changes



class ExclusiveSampleSet(SequentialSampleSet):
    """
    ExclusiveSampleSet has no own samples, only inferred Sampled from the
    underlying SampleSets
    """

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

class MovePath(object):
    """
    Information Class to contain the actual route/path of moves taken in
    a complex (Path)Mover object

    The contains all information necessary to
    (a) turn the initial_sampleset into the final_sampleset
    (b) contains the result of all submovers used in the generation of
        the final SampleSet, this includes the order in which the submoves
        are called, if their result is the same as before (thus accepted)
        and what samples were apply to the current (intermediate) SampleSet

    - The intermediate samplesets are not stored but used internally to track
    results.
    - Each submove has a SampleSet as result and samplesets are changed using
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

    SampleSets form a group like structure with our known concatenation of
    applying all samples in a row. But this requires that application of
    Samples is unique. Can we assure that? This requires to set an initial
    sample that is to be replaced by another or nothing. Adding just means to
    start from an emtpy sample. The randomness need to be stored in the sample.
    The we can really do

    setA * setB = setC

    setA * setB * setC = setA * (setB * setC) = (setA * setB) * setC

    setA * setB = []

    We need SampleSet to know the original SampleSet and a SampleSet
    that contain the moves, which can be nested. So A Move takes the full
    MoveSample and Returns the object with added Moves which are itself
    a SampleSet.

    SeqentialSampleSet(
        initial_set, [
        SampleSet(
            moveA1:y:initial_set + [sample1, sample2]
        ),
        PartialSampleSet(
        [
            moveB1:y:initial_set + [sample1, sample2] + [sample3, sample4],
            moveB2:n:initial_set + [sample1, sample2] + [sample3, sample4] + [sample5, sample6]
        ]):y:initial_set + [sample1, sample2] + [sample3, sample4],
        ExcluiveSampleSet([
            moveC1:y:initial_set + [sample1, sample2] + [sample3, sample4] + [sample7, sample8, sample9],
            moveC2:n:initial_set + [sample1, sample2] + [sample3, sample4] + [sample7, sample8, sample9] + [sample10, sample11]
        ]):n:[],
        moveA2:y:initial_set + [sample1, sample2] + [sample3, sample4] + [sample12]
    ]):y:initial_set + [sample1, sample2] + [sample3, sample4] + [sample12]


    Attributes
    ----------
    initial_sampleset : SampleSet
        the initial SampleSet where the path originates
    final_sampleset : SampleSet
        the final SampleSet where the path ends
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
        new SampleSet from the old one. Overwritten samples are removed.

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
