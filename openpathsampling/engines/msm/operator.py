import numpy as np


class O(object):
    """
    Define that a trajectory of a length in `lengths` was observed in the
    `states` distribution

    Example: If I say lengths: slice(None,None), states = [stateA: 100%] it
    means you assume that you were with 100% in stateA for an arbitrary number
    of steps
    """

    def __sub__(self, other):
        """
        Full Difference. Self & other = other

        Parameters
        ----------
        other

        Returns
        -------

        """

        return Diff(self, other)

    def __add__(self, other):
        """
        Disjoint union observations. Self & other = empty

        Parameters
        ----------
        other

        Returns
        -------

        """
        return Union([self, other])

    def __mul__(self, other):
        """
        Chain observations
        Parameters
        ----------
        other

        Returns
        -------

        """
        if isinstance(other, Chain):
            return Chain([self] + other.oos)
        else:
            return Chain([self, other])

    def as_matrix(self, model):
        return None

    def __call__(self, model):
        np.dot(np.dot(model.init, self.as_matrix(model)), model.one)
        return None


class State(O):
    """
    Define that a trajectory of a length in `lengths` was observed in the
    `states` distribution

    Example: If I say lengths: slice(None,None), states = [stateA: 100%] it
    means you assume that you were with 100% in stateA for an arbitrary number
    of steps
    """

    def __init__(self, length, ins, outs, inverted=False):
        self.block_ins = ins
        self.block_outs = outs
        self.length = length
        self.inverted = inverted

    def __str__(self):
        return 'S[%s,%s,%s,%s]' % (
            str(self.block_ins),
            str(self.block_outs),
            str(self.length),
            str(self.inverted)
        )

    def __add__(self, other):
        """
        Disjoint union observations. Self & other = empty

        Parameters
        ----------
        other

        Returns
        -------

        """
        if isinstance(other, State):
            # if self.block == other.block:
            #     if not self.length & other.length:
            #         return State(self.block, self.length | other.length)
            #
            # if self.length == other.length:
            #     if not self.block & other.block:
            #         return State(self.block | other.block, self.length)
            #

            ls = self.length & - other.length
            lo = other.length & - self.length
            lb = self.length & other.length

            return State([
                self,
                other,
                State(
                    self.block & other.block, self.length & other.length, True)
            ])

        return super(State, self).__add__(other)

    def __and__(self, other):
        if isinstance(other, State):
            if self.inverted and other.inverted:
                return State(
                    self.block | other.block, self.length | other.length, True)
            elif self.inverted:
                return other & self
            elif other.inverted:
                return Diff(
                    self,
                    State(other.block - self.block, self.length & other.length))
            else:
                return State(
                    self.block & other.block, self.length & other.length)
        else:
            raise ValueError('& and | only work for States')

    def __or__(self, other):
        if isinstance(other, State):
            return State(self.block | other.block, self.length | other.length)
        else:
            raise ValueError('& and | only work for States')

    def as_matrix(self, model):
        single_step = np.sum(model.basis[self.block.states])

        return self.length.matrix_mult(single_step)

    def __invert__(self):
        return State(self.block, self.length, ~ self.inverted)


class Chain(O):
    def __init__(self, oos):
        super(Chain, self).__init__()
        self.oos = oos

    def as_matrix(self, model):
        return reduce(np.dot, (o.as_matrix() for o in self.oos))


class Union(O):
    def __init__(self, oos):
        super(Union, self).__init__()
        self.oos = oos

    def as_matrix(self, model):
        return reduce(np.sum, (o.as_matrix() for o in self.oos))


class Diff(O):
    def __init__(self, o1, o2):
        super(Diff, self).__init__()
        self.o1 = o1
        self.o2 = o2

    def as_matrix(self, model):
        return self.o1.as_matrix(model) - self.o2.as_matrix(model)


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
        this.basis = np.array(
            np.diagonal(a) for a in msm
        )

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
#   -> Seq([ Length(1) & AllInXVolume(), AllOutAB & PartOutI,
# Length(1) & AllInXVolume(AB) ])
#   -> Seq([ State(A, 1), State(notAB, None) & not State(I, None), State(AB, 1])
#   -> Seq([ (A, 1), (notAB, None) & not (I, None), (AB, 1) ])
#   -> Seq([ (A, 1), (notAB, None) - (I & not AB, None), (AB, 1) ])
#   -> Seq([ (A, 1), (notAB, None) - (I, None), (AB, 1) ])


# MINUS = Sequential([Single(All(A)), AllOutA & PartOutI, Single(all(A or B))])
#   -> Seq([ Length(1) & AllInXVolume(), AllOutAB & PartOutI,
# Length(1) & AllInXVolume(AB) ])
#   -> Seq([ State(A, 1), State(notAB, None) & not State(I, None), State(AB, 1])
#   -> Seq([ (A, 1), (notAB, None) & not (I, None), (AB, 1) ])
#   -> Seq([ (A, 1), (notAB, None) - (I & not AB, None), (AB, 1) ])
#   -> Seq([ (A, 1), (notAB, None) - (I, None), (AB, 1) ])

# State(X, None) & not State(Y, None) -> State(X, None) - State(Y - X, None)
# State(X, L) & not State(Y, K) -> State(X, L) - State(Y - X, L & K)

# State(X, L) | not State(Y, K) -> State(X, L) + 1.0 - State(X | Y, K | L)


# State(X, None) & State(Y, None) -> State(Y & X, None)
# State(X, None) | State(Y, None) -> State(Y | X, None)

# State(X, L) & State(Y, K) -> State(Y & X, L & K)
# State(X, L) | State(Y, K) -> State(Y | X, L | K)
