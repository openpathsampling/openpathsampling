__author__ = 'jan-hendrikprinz'

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
