import numpy as np


class MSM(object):
    def __init__(self, matrix):
        self.matrix = matrix



class Lengths(object):
    # def __init__(self, inverted=False):
    #     self.inverted = inverted

    def __and__(self, other):
        return self

    def __or__(self, other):
        return self

    def __lt__(self, other):
        return all(l in self for l in other.lengths)

    def __gt__(self, other):
        return all(l in other for l in self.lengths)

    def __eq__(self, other):
        return self > other > self

    def __contains__(self, item):
        return False

    def __nonzero__(self):
        return self > EmptyLengths()

    @property
    def lengths(self):
        return []


class SetLengths(Lengths):
    def __init__(self, length_set):
        super(SetLengths, self).__init__()
        self.length_set = length_set

    def __nonzero__(self):
        return bool(self.length_set)

    def __and__(self, other):
        if other < self:
            return other
        elif other > self:
            return self

        if isinstance(other, IntLengths):
            return EmptyLengths()

        elif isinstance(other, SliceLengths):
            if other.length_slice.stop is not None:
                return SetLengths(
                    self.length_set &
                    set(range(
                        *other.length_slice.indices(other.length_slice.stop))))
            else:
                return SetLengths(
                    self.length_set &
                    set(range(
                        *other.length_slice.indices(max(self.length_set)))))
        elif isinstance(other, SetLengths):
            return SetLengths(
                self.length_set & other.length_set
            )
        elif isinstance(other, MultiLengths):
            return other & self

    def __or__(self, other):
        if other < self:
            return self
        elif other > self:
            return other

        if isinstance(other, IntLengths):
            return SetLengths({other.length} | self.length_set)

        elif isinstance(other, SliceLengths):
            if other.length_slice.stop is not None:
                return MultiLengths([
                    SetLengths(
                        self.length_set -
                        set(range(
                            *other.length_slice.indices(other.length_slice.stop)))),
                    other
                ])
            else:
                if min(self.length_set) >= other.length_slice.start:
                    return other
                else:
                    return MultiLengths([
                        SetLengths(
                            self.length_set -
                            set(range(
                                *other.length_slice.indices(max(self.length_set))))),
                        other
                    ])
        elif isinstance(other, SetLengths):
            return SetLengths(
                self.length_set | other.length_set
            )
        elif isinstance(other, MultiLengths):
            return other | self

    def __contains__(self, item):
        if isinstance(item, int):
            return item in self.length_set

        elif isinstance(item, slice):
            if item.stop is not None:
                return set(range(*item.indices(item.stop))) <= self.length_set
            else:
                return False
        elif isinstance(item, set):
            return item <= self.length_set

    @property
    def lengths(self):
        return [self.length_set]


class MultiLengths(Lengths):
    def __init__(self, length_list):
        super(MultiLengths, self).__init__()
        self.length_list = [l for l in length_list if l]

    def __contains__(self, item):
        for l in self.length_list:
            if item in l:
                return True

        return False

    @property
    def lengths(self):
        return self.length_list

    def __nonzero__(self):
        return any(self.length_list)


    def __and__(self, other):
        if other < self:
            return other
        elif other > self:
            return self

        return MultiLengths([
            l & self for l in other.length_list
        ])

    def __or__(self, other):
        if other < self:
            return self
        elif other > self:
            return other

        return MultiLengths([
            l | self for l in other.length_list
        ])


class IntLengths(Lengths):
    def __init__(self, length):
        super(IntLengths, self).__init__()
        self.length = length

    def __nonzero__(self):
        return True

    def __contains__(self, item):
        if isinstance(item, int):
            return item == self.length
        elif isinstance(item, slice):
            if self.length == 0:
                if item.start is None or item.start == 0 and item.stop == 1:
                    return True
            else:
                if item.start == self.length and item.stop == self.length + 1:
                    return True

            return False
        elif isinstance(item, set):
            return {self.length} == item

    def __and__(self, other):
        if other > self:
            return self

        # if the single state is not smaller than other it must be empty
        return EmptyLengths()

    def __or__(self, other):
        if other < self:
            return self

        elif other > self:
            return other

        if isinstance(other, IntLengths):
            return SetLengths({other.length, self.length})

        elif isinstance(other, SliceLengths):
            if self.length < other.length_slice.start or (other.length_slice.stop is not None and self.length > other.length_slice.stop - 1):
                return MultiLengths([
                    self,
                    other
                ])
            else:
                return other
        elif isinstance(other, SetLengths):
            return SetLengths(
                {self.length} | other.length_set
            )
        elif isinstance(other, MultiLengths):
            return other | self

    @property
    def lengths(self):
        return [self.length]


class SliceLengths(Lengths):

    def __nonzero__(self):
        if self.length_slice.stop is None:
            return True
        elif self.length_slice.start is None:
            return self.length_slice.stop is not 0
        else:
            return self.length_slice.start < self.length_slice.stop

    def __init__(self, length_slice):
        super(SliceLengths, self).__init__()
        self.length_slice = length_slice

    def __contains__(self, item):
        if isinstance(item, int):
            if (self.length_slice.start is None or self.length_slice.start <= item) and \
                (self.length_slice.stop is None or self.length_slice.stop > item):
                    return True

            return False
        elif isinstance(item, slice):
            if item.start is None and item.stop is None:
                return self.length_slice.start is None and self.length_slice.stop is None
            elif item.start is None:
                return self.length_slice.start is None and self.length_slice.stop >= item.stop
            elif item.stop is None:
                return self.length_slice.stop is None and self.length_slice.start <= item.start
            else:
                return self.length_slice.start <= item.start and self.length_slice.stop >= item.stop

        elif isinstance(item, set):
            return self.length_slice.start <= min(item) and self.length_slice.stop > max(item)

    def __and__(self, other):
        if other < self:
            return other
        elif other > self:
            return self

        if isinstance(other, IntLengths):
            return other & self

        elif isinstance(other, SliceLengths):
            if other.length_slice.stop is None:
                stop = self.length_slice.stop
            elif self.length_slice.stop is None:
                stop = other.length_slice.stop
            else:
                stop = min(self.length_slice.stop, other.length_slice.stop)

            if other.length_slice.start is None:
                start = self.length_slice.start
            elif self.length_slice.start is None:
                start = other.length_slice.start
            else:
                start = min(self.length_slice.start, other.length_slice.start)

            if start is not None and stop is not None:
                stop = max(start, stop)

            return SliceLengths(slice(start, stop))
        elif isinstance(other, SetLengths):
            return other & self
        elif isinstance(other, MultiLengths):
            return other & self

    def __or__(self, other):
        if other < self:
            return self
        elif other > self:
            return other

        if isinstance(other, IntLengths):
            return other | self

        elif isinstance(other, SliceLengths):
            if other.length_slice.stop is None:
                stop = None
                left = self.length_slice.stop
            elif self.length_slice.stop is None:
                stop = None
                left = other.length_slice.stop
            else:
                stop = max(self.length_slice.stop, other.length_slice.stop)
                left = min(self.length_slice.stop, other.length_slice.stop)

            if other.length_slice.start is None:
                start = None
                right = self.length_slice.start
            elif self.length_slice.start is None:
                start = None
                right = other.length_slice.start
            else:
                start = max(self.length_slice.start, other.length_slice.start)
                right = min(self.length_slice.start, other.length_slice.start)

            if left is not None and right is not None and left < right:
                return MultiLengths([
                    self, other
                ])
            else:
                return SliceLengths(slice(start, stop))
        elif isinstance(other, SetLengths):
            return other | self
        elif isinstance(other, MultiLengths):
            return other | self


    @property
    def lengths(self):
        return [self.length_slice]

class EmptyLengths(Lengths):
    def __init__(self):
        super(EmptyLengths, self).__init__(None)

    def __nonzero__(self):
        return False


    def __and__(self, other):
        return self

    def __or__(self, other):
        return other

    def __contains__(self, item):
        return False

    @property
    def lengths(self):
        return []


class AllLengths(SliceLengths):
    def __init__(self):
        super(AllLengths, self).__init__(None, None)

    def __nonzero__(self):
        return True

    def __and__(self, other):
        return other

    def __or__(self, other):
        return self

    def __contains__(self, item):
        return True


class OneLength(IntLengths):
    def __init__(self):
        super(OneLength, self).__init__(1)

    def __and__(self, other):
        if not isinstance(other, EmptyLengths)
            return

        return NotImplemented

    def __or__(self, other):
        return other

class FixedLength(IntLengths):
    def __init__(self, length):
        super(FixedLength, self).__init__(length)

class RangeLength(SliceLengths):
    def __init__(self, min_length, max_length):
        super(RangeLength, self).__init__(slice(min_length, max_length))

class UnionLength(MultiLengths):
    def __init__(self, lengths):
        super(UnionLength, self).__init__(sum((l.lengths for l in lengths), []))


class O(object):
    """
    Define that a trajectory of a length in `lengths` was observed in the
    `states` distribution

    Example: If I say lengths: slice(None,None), states = [stateA: 100%] it
    means you assume that you were with 100% in stateA for an arbitrary number
    of steps
    """

    def __init__(self, block, length):
        self.block = block
        self.length = length

    def __sub__(self, other):
        """
        Full Difference. Self & other = other

        Parameters
        ----------
        other

        Returns
        -------

        """


    def __add__(self, other):
        """
        Disjoint union observations. Self & other = empty

        Parameters
        ----------
        other

        Returns
        -------

        """
        return Joined()

    def __mul__(self, other):
        """
        Chain observations
        Parameters
        ----------
        other

        Returns
        -------

        """
        return Chain([self, other])

    def __and__(self, other):
        return O(self.block & other.block, self.length & other.length)

    def __or__(self, other):
        return O(self.block | other.block, self.length | other.length)

    def as_matrix(self, model):
        return None

    def __call__(self, model):
        np.dot(np.dot(model.init, self.as_matrix(model), model.one))
        return None


class Single(O):
    def __init__(self, states):
        super(Single, self).__init__()
        self.states = states

    def as_matrix(self, model):
        """
        Use a linear combination to return tau_state for lengths

        Parameters
        ----------
        model

        Returns
        -------

        """
        return


class Multi(O):
    def __init__(self, obs, lengths):
        super(Multi, self).__init__()
        self.lengths = lengths


class Chain(O):
    def __init__(self, oos):
        super(Chain, self).__init__(None, None)
        self.oos = oos

    def as_matrix(self, model):
        """
        Product

        Parameters
        ----------
        model

        Returns
        -------

        """
        return None


class OOM(object):
    def __init__(self):
        self.one = 1
        self.init = 1
        self.basis = []
        self.states = []
        self.size = 1

    @staticmethod
    def from_msm(msm):
        this = OOM()
        _size = len(msm)
        this.states = msm.states
        this.init = np.ones(_size) / _size
        this.one = np.ones(_size)
        this.basis = [
            np.diagonal(a) for a in msm
        ]

        return this

    def __len__(self):
        return self.size

    def observe_state(self, state):
        pass

# oo = SequentialEnsemble([...]).oo
# oo(OOM.from_msm(msm)) = 0.2

# AllInX() -> O({state : 1.0}, len=None))
# Sequential() -> Chain([ens.oo for ens in self.ensembles]
# SingleFrame(AllInX()) = Length(1) & AllInXVolume()
#   -> State(All, len=1)) & StateX({A : 1.0}, len=None))
#   -> State({A : 1.0}, len=1)

# TIS = Sequential([Single(All(A)), AllOutA & PartOutI, Single(all(A or B))])
#   -> Seq([ Length(1) & AllInXVolume(), AllOutAB & PartOutI, Length(1) & AllInXVolume(AB) ])
#   -> Seq([ State(A, 1), State(notAB, None) & not State(I, None), State(AB, 1])
#   -> Seq([ (A, 1), (notAB, None) & not (I, None), (AB, 1) ])
#   -> Seq([ (A, 1), (notAB, None) - (I & not AB, None), (AB, 1) ])
#   -> Seq([ (A, 1), (notAB, None) - (I, None), (AB, 1) ])


# MINUS = Sequential([Single(All(A)), AllOutA & PartOutI, Single(all(A or B))])
#   -> Seq([ Length(1) & AllInXVolume(), AllOutAB & PartOutI, Length(1) & AllInXVolume(AB) ])
#   -> Seq([ State(A, 1), State(notAB, None) & not State(I, None), State(AB, 1])
#   -> Seq([ (A, 1), (notAB, None) & not (I, None), (AB, 1) ])
#   -> Seq([ (A, 1), (notAB, None) - (I & not AB, None), (AB, 1) ])
#   -> Seq([ (A, 1), (notAB, None) - (I, None), (AB, 1) ])

# State(X, None) & not State(Y, None) -> State(X, None) - State(Y - X, None)
# State(X, L) & not State(Y, K) -> State(X, L) - State(Y - X, K - L)

# State(X, L) | not State(Y, K) -> State(X, L) + 1.0 - State(X | Y, K | L)


# State(X, None) & State(Y, None) -> State(Y & X, None)
# State(X, None) | State(Y, None) -> State(Y | X, None)

# State(X, L) & State(Y, K) -> State(Y & X, L & K)
# State(X, L) | State(Y, K) -> State(Y | X, L | K)

