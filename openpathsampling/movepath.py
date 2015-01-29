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

    [moveA1:y, moveB1:y, moveB2:y, moveC1:y, moveC2:y, moveA2:y]
    [moveA1:y, moveB1:y, moveB2:y, moveC1:y, moveC2:n, moveA2:y]
    [moveA1:y, moveB1:y, moveB2:y, moveC1:n,           moveA2:y]
    [moveA1:y, moveB1:y, moveB2:n, moveC1:y, moveC2:y, moveA2:y]
    [moveA1:y, moveB1:y, moveB2:n, moveC1:y, moveC2:n, moveA2:y]
    [moveA1:y, moveB1:y, moveB2:n, moveC1:n,           moveA2:y]
    [moveA1:y, moveB1:n,           moveC1:y, moveC2:y, moveA2:y]
    [moveA1:y, moveB1:n,           moveC1:y, moveC2:n, moveA2:y]
    [moveA1:y, moveB1:n,           moveC1:n,           moveA2:y]

    S([moveA1:y, P([moveB1:y, moveB2:y]):y, E([moveC1:y, moveC2:y]):y, moveA2:y]):y
    S([moveA1:y, P([moveB1:y, moveB2:y]):y, E([moveC1:y, moveC2:n]):n, moveA2:y]):y
    S([moveA1:y, P([moveB1:y, moveB2:y]):y, E([moveC1:n]):n,           moveA2:y]):y
    S([moveA1:y, P([moveB1:y, moveB2:n]):y, E([moveC1:y, moveC2:y]):y, moveA2:y]):y
    S([moveA1:y, P([moveB1:y, moveB2:n]):y, E([moveC1:y, moveC2:n]):n, moveA2:y]):y
    S([moveA1:y, P([moveB1:y, moveB2:n]):y, E([moveC1:n]):n,           moveA2:y]):y
    S([moveA1:y, P([moveB1:n]):y,           E([moveC1:y, moveC2:y]):y, moveA2:y]):y
    S([moveA1:y, P([moveB1:n]):y,           E([moveC1:y, moveC2:n]):n, moveA2:y]):y
    S([moveA1:y, P([moveB1:n]):y,           E([moveC1:n]):n,           moveA2:y]):y

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
