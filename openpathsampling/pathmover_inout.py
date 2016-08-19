from collections import Counter
from itertools import product


# ------------------------------------------------------------------------------
#   Replica In-Out-Logic
# ------------------------------------------------------------------------------

# The following classes are used to inspect the effects of a PathMover
# on the input SampleSet w.r.t. the output SampleSet
# As an example: A ReplicaExchangePathMover will switch one replica from
# ens1 and exchange it with one replica from ens2. We could express this
# as ens1 -> ens2 and ens2 -> ens1. Every mover can be expressed in this
# way with possible multiple occurrences of a move of a replica

# at the time this cannot handle movers that pick the used ensembles conditioned
# on the actual sample_set that means `FirstAllowedMover`, `LastAllowedMover`
# and `RandomAllowedChoiceMover` are not inspected properly. It will still
# give the potential list of all possible InOuts but using conditions this
# can be smaller. The general change of this to graph-based analysis will
# be done in 2.0

# Currently this feature is only used for SRTIS and for some kinds of
# _bootstrapping_ / generation of initial samples


class ReplicaStateSet(set):
    """
    Represents a set of possible state of replicas

    See Also
    --------
    `ReplicaState`, `InOut`, `InOutSet`


    """

    @staticmethod
    def from_sampleset(sample_set):
        """
        Construct a set of a single state from a `SampleSet`

        Parameters
        ----------
        sample_set : :obj:`openpathsampling.SampleSet`
            The sampleset turned into a single set replica state

        Returns
        -------
        :obj:`ReplicaStateSet`
            the constructed set of replica states

        """
        return ReplicaStateSet({ReplicaState.from_sampleset(sample_set)})

    @staticmethod
    def from_ensembles(ensembles):
        """
        Construct a set of a single state from a list of ensembles

        Parameters
        ----------
        ensembles : iterable of :obj:`openpathsampling.Ensemble`
            The ensembles turned into a single set replica state

        Returns
        -------
        :obj:`ReplicaStateSet`
            the constructed set of replica states

        """
        return ReplicaStateSet({ReplicaState.from_ensembles(ensembles)})

    @staticmethod
    def from_ensembles_dict(ensembles_dict):
        """
        Construct a set of a single state from a dictionary of ensembles to ints

        Parameters
        ----------
        ensembles_dict : dict of :obj:`openpathsampling.Ensemble`: int
            The dict representing the number of times an ensemble is in the
             replica state

        Returns
        -------
        :obj:`ReplicaStateSet`
            the constructed set of replica states

        """
        return ReplicaStateSet(
            {ReplicaState.from_ensemble_dict(ensembles_dict)})

    def _reduce(self, func):
        return reduce(func, map(lambda x: Counter(dict(x)), self), Counter())

    def reduce_min(self):
        return self._reduce(lambda x, y: x & y)

    def reduce_max(self):
        return self._reduce(lambda x, y: x & y)


class ReplicaState(frozenset):
    """
    Represents a set of samples: how many samples per ensembles

    This object is represented bya frozenset of tuples to make it hashable.
    Technically it could be represented by a :class:`collections.Counter`

    lesser than and greater than are implemented and work as they would for
    a Counter. So lesser or equal means that all ensembles present in the
    "smaller" state are also present in the "larger" one and the
    multiplicity is smaller for all ensembles.

    This is useful to check if certain requirements are met. When `necessary`
    represent the minimal necessary number of samples per ensemble and `current`
    is the current state of the sample_set then `necessary <= current` checks if
    the requirements are met

    Replica states allow comparison with inclusion using `>` and `<`. So, if
    all ensembles from A are also present in B and all of the multiplicities of
    A are smaller than that of B then A < B

    See Also
    --------
    `ReplicaStateSet`, `InOut`, `InOutSet`

    """

    @staticmethod
    def from_sampleset(sample_set):
        """
        Construct a `ReplicaState` from a sampleset

        Parameters
        ----------
        sample_set : `openpathsampling.SampleSet`
            the sampleset to be condensed into a ReplicaState

        Returns
        -------
        `ReplicaState`
            the replicastate representing the multiplicity in ensembles
            present in the sampleset
        """
        d = {}
        for sample in sample_set:
            d[sample.ensemble] = d.get(sample.ensemble, 0) + 1

        return ReplicaState(d.items())

    @staticmethod
    def from_ensembles(ensembles):
        """
        Construct a `ReplicaState` from a list of ensembles

        Parameters
        ----------
        ensembles :list of  `openpathsampling.Ensemble`
            the list of ensembles to be turned into a ReplicaState

        Returns
        -------
        `ReplicaState`
            the replicastate representing the multiplicity in ensembles
            present in the sampleset. In this case each ensemble is used
            with multiplicity one.
        """

        d = {}
        for ens in ensembles:
            d[ens] = 1

        return ReplicaState(d.items())

    @staticmethod
    def from_ensemble_dict(ensemble_dict):
        """
        Construct a `ReplicaState` from a ensemble dictionary

        Parameters
        ----------
        ensemble_dict : dict(`openpathsampling.Ensemble`: int)
            the dictionary turned into a replica state. keys are the
            ensembles used and the value is the multiplicity

        Returns
        -------
        `ReplicaState`
            the replicastate representing the multiplicity in ensembles
            present in the sampleset
        """

        return ReplicaState(ensemble_dict.items())

    def __str__(self):
        ensemble_list = sorted([s for s in self], key=lambda x: hex(id(x[0])))

        s = []
        for ens in ensemble_list:
            s += ["{:>30} ({:>11}) : {:>3}".format(
                ens[0].name, hex(id(ens[0])), ens[1])]

        return '\n'.join(s)

    def filter(self, ensembles):
        """
        Filter a replica state by a list of ensembles

        Parameters
        ----------
        ensembles : list of `openpathsampling.Ensembles`
            the list of ensembles which represent the filter. Only ensembles
            that are als in `ensembles` we be kept with their respective
            multipliity

        Returns
        -------
        `ReplicaState`
            the reduced filtered replica state
        """
        return ReplicaState({s for s in self if s[0] in ensembles})

    def __gt__(self, other):
        return Counter(dict(self)) > Counter(dict(self))

    def __lt__(self, other):
        return Counter(dict(self)) < Counter(dict(self))


class InOutSet(set):
    """
    Represents a set of possible in-out relations

    See Also
    --------
    `ReplicaState`, `ReplicaStateSet`, `InOut`

    """

    def __add__(self, other):
        if other is None or len(other) == 0:
            return self
        elif len(self) == 0:
            return other
        else:
            return InOutSet(set.union(*[
                in1 * in2
                for in1 in self for in2 in other
            ]))

    def __radd__(self, other):
        if other is None or len(other) == 0:
            return self
        elif len(self) == 0:
            return other
        else:
            return InOutSet(set.union(*[
                in1 * in2
                for in1 in other for in2 in self
            ]))

    @property
    def ins_minimal(self):
        c = Counter()
        for s in self:
            if s.essential:
                c |= s.ins

        return c

    @property
    def ins(self):
        """
        The maximally needed replica state for input

        Maximally means, that larger input will not change the behaviour
        anymore or cause different behaviour

        A mover might be called "simple" if the minimal and maximally required
        replica state is the same. We could additionally require that the
        multiplicity per ensemble is one.

        Returns
        -------
        :obj:'collections.Counter`
            a Counter object representing the maximal replica state used
            for input

        Notes
        -----
        A counter can be turned into a ReplicaState by
        `ReplicaState(dict(counter).items())`.

        """
        c = Counter()
        for s in self:
            c |= s.ins

        return c

    @property
    def outs(self):
        """
        The maximal set of ensembles
        Returns
        -------

        """
        c = Counter()
        for s in self:
            c |= s.outs

        return c

    @property
    def outs_minimal(self):
        """

        Returns
        -------

        """
        c = Counter()
        for s in self:
            if s.essential:
                c |= s.outs

        return c

    @property
    def is_constant(self):
        """
        Check whether the move will keep the number of samples per ensemble

        Returns
        -------
        bool

        """
        return all([s.ins == s.outs for s in self])

    def filter(self, ensembles):
        """
        Return InOutSet with relations within a given set of ensembles

        Parameters
        ----------
        ensembles : iterable of `openpathsampling.Ensemble`

        Returns
        -------
        `InOutSet`
            the reduced in-out-relation table
        """
        return InOutSet({
            s.filter(ensembles) for s in self if set(s.ins) <= set(ensembles)
        })

    def move(self, replica_states):
        """
        Move a set of replica states and return a set of possible outcomes

        Parameters
        ----------
        replica_states : `ReplicaStateSet`

        Returns
        -------
        `ReplicaStateSet`
            the set of possible replica states being produced by this
            in-out-relation
        """
        ret = set()

        for replica_state in replica_states:
            c = Counter(replica_state)

            # Only move a replica state if it meets the minimal requirements
            # are met.
            if c >= self.ins_minimal:
                ret.update({s.move(replica_state) for s in self if s.ins <= c})

        return ReplicaStateSet(ret)

    def filter_max_length(self, length):
        """
        Return InOutSet limited by maximal number of samples per ensemble

        Parameters
        ----------
        length : int
            the number of ensembles

        Returns
        -------
        `InOutSet`
        the reduced in-out-relation table
        """
        return InOutSet({s for s in self if max(s.ins.values()) <= length})

    def mixing_matrix(self, ensembles):
        """
        Return a matrix of bool if a sample can move between ensembles

        Returns
        -------

        """
        res = {
            ens1: {
                ens2: set()
                for ens2 in ensembles
            } for ens1 in ensembles
        }
        for s in self:
            for (e1, e2, fix), v in s:
                res[e1][e2].add(fix)

        return res


class InOut(frozenset):
    """
    Represent the change in occupied ensembles during a move

    A mover changes a sampleset and will replace samples, move
    them between ensembles or create new ones. This change will be
    represented by this object.

    Assume that we replace a sample in ensemble A and move a second
    sample from A to B. This wil be represented by

        { (A, A, 1),
          (A, B, 1) }

    This class makes handling these objects easier. Chaining these changes
    will result in new possible changes. Since chaining can result in multiple
    possible relations depending on the occupation state of the ensemble
    chaingin will return a set of in-out-relations.
    """

    _n_max_samples_per_ensemble = 1
    _use_move_type = 0

    def __new__(cls, *args):
        return frozenset.__new__(cls, args[0])

    def __init__(self, relations=None, essential=None):
        # note that frozenset uses __new__ to input data. The next line is
        # only for making sure all the rest is set correctly.
        # __init__ will NOT set the content!
        frozenset.__init__(self, relations)

        if essential is None:
            essential = True

        self.essential = essential

    @property
    def ensembles(self):
        """
        Return a list of all appearing ensembles in the relations

        Returns
        -------

        """
        outs = set([s[1] for s in self])
        ins = set(r[0] for r in self)
        ens = ins | outs

        return ens

    @property
    def ins(self):
        """
        Return input ensembles and their multiplicity

        Returns
        -------
        :obj:'collections.Counter`
            a Counter object representing the input requirements

        """
        d = Counter()
        for s, v in self:
            d[s[0]] += v

        return d

    @property
    def outs(self):
        """
        Returns
        -------
        :obj:'collections.Counter`
            a Counter object representing the input requirements

        """

        d = Counter()
        for s, v in self:
            d[s[1]] += v

        return d

    def filter(self, ensembles):
        """
        Remove all relations in using a set of ensembles

        Parameters
        ----------
        ensembles : iterable of `openpathsampling.Ensemble`
            the set of ensembles

        Returns
        -------

        """
        return InOut(
            [s for s in self if s[0][0] in ensembles and s[0][1] in ensembles])

    # def in_ensembles(self, ensembles):
    #     return not bool(self.ensembles - set(ensembles))

    def __mul__(self, other):

        mat1 = dict(self)
        mat2 = dict(other)

        outs = set([s[1] for s in mat1])
        ins = set(r[0] for r in mat2)
        ens = ins | outs

        parts = []

        # do this for all inner ensembles
        for e in ens:
            # now we have `froms -> e -> tos` as all possibilities with
            # one connection and we need to form all possible pair combinations
            # with zero, one, ... pairs
            froms = sum([[(s[0], s[2])] * mat1[s] for s in mat1 if s[1] is e], [])
            tos = sum([[(s[1], s[2])] * mat2[s] for s in mat2 if s[0] is e], [])

            parts.append(self._fromto(froms, e, tos))

        if self.essential and other.essential:
            outs = InOutSet(
                map(lambda x: InOut(
                    sum(zip(*x)[0], Counter()).items(), all(zip(*x)[1])),
                    product(*parts)))
        else:
            outs = InOutSet(
                map(lambda x: InOut(
                    sum(zip(*x)[0], Counter()).items(), False),
                    product(*parts)))

        if InOut._n_max_samples_per_ensemble > 0:
            return outs.filter_max_length(InOut._n_max_samples_per_ensemble)
        else:
            return outs

    def _fromto(self, froms, e, tos):
        frees = [(f[0], e, f[1]) for f in froms] + [(e, t[0], t[1]) for t in tos]
        pairs = list(product(range(len(froms)), range(len(tos))))
        if len(pairs) == 0:
            inouts = [(Counter(frees), True)]
        else:
            inouts = list()
            inouts.append((Counter(frees), False))

            for i1, i2 in pairs:
                fix = (froms[i1][0], tos[i2][0], froms[i1][1] * tos[i2][1] * InOut._use_move_type)
                r_froms = froms[:i1] + froms[i1 + 1:]
                r_tos = tos[:i2] + tos[i2 + 1:]
                for rest in self._fromto(r_froms, e, r_tos):
                    if fix[2] < 1 or fix[0] is not fix[1]:
                        rest[0][fix] += 1

                    if rest[0][fix] > 0:
                        inouts.append(rest)

        return inouts

    def move(self, replica_state):
        """
        Move a replica state by the relations in this object

        Parameters
        ----------
        replica_state

        Returns
        -------

        """
        d = Counter(dict(replica_state))
        d = d - self.ins + self.outs

        return ReplicaState(d.items())