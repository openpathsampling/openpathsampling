import itertools
import collections
import abc

import openpathsampling as paths
from openpathsampling.netcdfplus import StorableNamedObject

from functools import reduce  # not built in for py3

LevelLabels = collections.namedtuple(
    "LevelLabels",
    ["SIGNATURE", "MOVER", "GROUP", "SUPERGROUP", "GLOBAL"]
)


def most_common_value(ll):
    """
    Calculates the most common value and its count.
    """
    return collections.Counter(ll).most_common(1)[0]


class StrategyLevels(LevelLabels):
    """
    Custom version of a namedtuple to handle aspects of the `level`
    attribute of MoveStategy.
    """
    def level_type(self, lev):
        """
        Determines which defined level the value `lev` is closest to. If
        the answer is not unique, returns `None`.
        """
        if lev < 0 or lev > 100:
            return None
        levels = [self.SIGNATURE, self.MOVER, self.GROUP, self.SUPERGROUP,
                  self.GLOBAL]
        distances = [abs(lev - v) for v in levels]
        mindist = min(distances)
        indices = [i for i in range(len(distances))
                   if distances[i] == mindist]
        if len(indices) > 1:
            return None
        else:
            return levels[indices[0]]


# possible rename to make available as paths.strategy_levels?
levels = StrategyLevels(
    SIGNATURE=10,
    MOVER=30,
    GROUP=50,
    SUPERGROUP=70,
    GLOBAL=90
)


class MoveStrategy(StorableNamedObject):
    """
    Each MoveStrategy describes one aspect of the approach to the overall
    MoveScheme. Within path sampling, there's a near infinity of reasonable
    move schemes to be used; we use MoveStrategy to simplify mixing and
    matching of different approaches.

    Parameters
    ----------
    ensembles : list of list of Ensemble, list of Ensemble, Ensemble,  or None
        The ensembles used by this strategy
    group : string or None
        The group this strategy is associated with (if any).
    replace : bool
        Whether this strategy should replace existing movers. See also
        `MoveStrategy.set_replace` and `MoveScheme.apply_strategy`.

    Attributes
    ----------
    replace_signatures : bool or None
        Whether this strategy should replace at the signature level.
    replace_movers : bool or None
        Whether this strategy should replace at the mover level.
    level
    """
    _level = -1

    __metaclass__ = abc.ABCMeta

    def __init__(self, ensembles, group, replace):
        super(MoveStrategy, self).__init__()
        self.ensembles = ensembles
        self.group = group
        self.replace = replace
        self.replace_signatures = None
        self.replace_movers = None
        self.set_replace(replace)
        self.from_group = group

    @property
    def level(self):
        """
        The level of this strategy.

        Levels are numeric, but roughly correspond to levels in the default
        move tree. This way, we build the tree from bottom up.
        """
        return self._level

    @level.setter
    def level(self, value):
        self._level = value
        self.set_replace(self.replace)  # behavior depends on level

    def set_replace(self, replace):
        """Sets values for replace_signatures and replace_movers."""
        self.replace_signatures = False
        self.replace_movers = False
        self.replace_group = False
        level_type = levels.level_type(self.level)
        if level_type == levels.SIGNATURE:
            self.replace_signatures = replace
        elif level_type == levels.MOVER:
            self.replace_movers = replace
        elif level_type == levels.SUPERGROUP:
            self.replace_group = replace

    def get_ensembles(self, scheme, ensembles):
        """
        Regularizes ensemble input to list of list.

        Input None returns a list of lists for ensembles in each sampling
        transition in the network. A list of ensembles is wrapped in a list.

        Parameters
        ----------
        ensembles : None, list of Ensembles, or list of list of Ensembles
            Input ensembles.

        Returns
        -------
        list of list of Ensembles
            Regularized output.

        Note
        ----
            List-of-list notation is used, as it is the most generic, and
            likely to be useful for many types of strategies. If desired,
            a strategy can always flatten this after the fact.
        """
        if ensembles is None:
            res_ensembles = []
            for t in scheme.network.sampling_transitions:
                res_ensembles.append(t.ensembles)
        else:
            # takes a list and makes it into list-of-lists
            res_ensembles = []
            elem_group = []
            try:
                ens_iter = iter(ensembles)
            except TypeError:
                ens_iter = iter([ensembles])
            for elem in ens_iter:
                try:
                    append_group = list(elem)
                except TypeError:
                    elem_group.append(elem)
                else:
                    if len(elem_group) > 0:
                        res_ensembles.append(elem_group)
                        elem_group = []
                    res_ensembles.append(append_group)
            if len(elem_group) > 0:
                res_ensembles.append(elem_group)

        return res_ensembles

    def get_per_mover_ensembles(self, scheme):
        """Get ensembles for each mover to be created.

        Every mover has specific initialization ensembles. This method
        figures out which ones to use for each mover to be created, based on
        the input scheme (and the network in it) and the ``self.ensembles``
        stored in this strategy.

        Parameters
        ----------
        scheme : :class:`.MoveScheme`
            move scheme with network to be used

        Returns
        -------
        list of list of :class:`.Ensemble` :
            the ensembles for each mover to be created: outer list is over
            movers; inner list is over ensembles for the given mover
        """
        raise NotImplementedError

    def get_parameters(self, scheme, list_parameters=None,
                       nonlist_parameters=None):
        """
        Gather initialization parameters for each mover to be created.

        Parameters
        ----------
        scheme : :class:`.MoveScheme`
            move scheme from which we determine the ensembles required by
            each mover
        list_parameters : list
            each item in the list can be either (1) a non-iterable item,
            in which case the same item will be used for all movers; or (2)
            a list of items, which correspond one-to-one with the movers to
            be created. Use this for parameters that might vary depending on
            the mover.
        nonlist_parameters : list
            each item in this list will be used for all movers. Use this for
            parameters that are the same for all movers (e.g., same engine
            for all shooting movers)

        Returns
        -------
        list of list :
            list containing the list of parameters for each mover; the
            specific strategy may need to unpack in substructure in its
            make_movers method. Order is: ensembles (in order returned by
            :meth:`.get_per_mover_ensembles`), list_parameters (in order
            given as input), nonlist_parameters (in order given as input)
        """
        ensemble_list = self.get_per_mover_ensembles(scheme)
        n_movers = len(ensemble_list)

        if list_parameters is None:
            list_parameters = []
        if nonlist_parameters is None:
            nonlist_parameters = []

        list_params = []
        for param in list_parameters:
            try:
                n_param = len(param)
            except TypeError:
                list_params.append([param] * n_movers)
            else:
                if n_param != n_movers:
                    raise RuntimeError(
                        "Error in move strategy parameters: found %d items "
                        + "for %d movers: %s", n_param, n_movers, str(param)
                    )
                else:
                    list_params.append(param)

        nonlist_params = [[param] * n_movers
                          for param in nonlist_parameters]

        all_params = list(zip(*ensemble_list)) + list_params + nonlist_params
        return list(zip(*all_params))

    @abc.abstractmethod
    def make_movers(self, scheme):
        """
        Makes the movers associated with this strategy.

        The exact behavior of this function differs somewhat depending on
        the `strategy.level`. In particular, this function can have
        side-effects on the `scheme`.

        For example, `GLOBAL`-level strategies must set
        `scheme.choice_probability`.

        Parameters
        ----------
        scheme : paths.MoveScheme
            the move scheme that this strategy will be used for

        Returns
        -------
        paths.PathMover or list of paths.PathMover
            the movers created by this part of the strategy
        """
        raise NotImplementedError # TODO: use ABCError when 302 is merged


class SingleEnsembleMoveStrategy(MoveStrategy):
    """Abstract class for movers with one (and same) input/output ensemble.
    """
    @staticmethod
    def signatures_to_init_ensembles(signatures):
        """
        Convert signatures into initialization ensembles.

        Note
        ----
            Subclasses may need to override this. The default behavior works
            for movers where the input and output have only one ensemble,
            and that ensemble is the same for both.

        Parameters
        ----------
        signatures : list of 2-tuple of tuple of :class:`.Ensemble`
            the move signatures from a path mover

        Returns
        -------
        list of :class:`.Ensemble`
            the ensembles to be used for initialization
        """
        init_ensembles = [s[0] for s in signatures]
        for e in init_ensembles: assert(len(e) == 1)  # sanity check
        return init_ensembles

    def get_init_ensembles(self, scheme):
        """
        Procedure for selecting the initialization ensembles.

        If this strategy has a specific set of ensembles, then those
        ensembles are used (after wrapping into list-of-list form).
        If the strategy does not have movers, but `self.from_group` already
        exists in the scheme, then we use the movers from the scheme to
        identify the signatures, and therefore the ensembles.
        If neither the strategy nor the scheme provide a list of ensembles,
        we use the default of all ensembles (initially grouped by
        transitions). If self.replace_signatures has been set to True, then
        we ignore the 

        Parameters
        ----------
        scheme : :class:`.MoveScheme`
            the move scheme this strategy will be applied to

        Returns
        -------
        list of list of :class:`.Ensemble`
            the ensembles to be used to initialize the movers
        """
        # Consider 3 cases:
        # if self.ensembles is not None, use those
        # elif scheme.movers[from_group], extract from that
        # elif self.ensembles is None, use all the ensembles (as now)
        try:
            movers = scheme.movers[self.from_group]
        except KeyError:
            movers = []

        # these are the conditions
        conditions = (self.ensembles is not None
                      or len(movers) == 0
                      or self.replace_signatures or self.replace_movers)

        if conditions:
            # covers cases 1 and 3 above
            init_ensembles = self.get_ensembles(scheme, self.ensembles)
        else:
            # covers case 2 above
            sigs = [m.ensemble_signature for m in movers]
            sig_ensembles = self.signatures_to_init_ensembles(sigs)
            init_ensembles = self.get_ensembles(scheme, sig_ensembles)

        return init_ensembles

    def get_per_mover_ensembles(self, scheme):
        ensemble_list = self.get_init_ensembles(scheme)
        ensembles = reduce(list.__add__, map(lambda x: list(x), ensemble_list))
        return [[ens] for ens in ensembles]



class OneWayShootingStrategy(SingleEnsembleMoveStrategy):
    """
    Strategy for OneWayShooting. Allows choice of shooting point selector.

    Parameters
    ----------
    selector : :class:`.ShootingPointSelector`
        method used to select shooting point
    ensembles : list of :class:`.Ensemble`
        ensembles for which this strategy applies; None gives default
        behavior
    engine : :class:`.DynamicsEngine`
        engine for the dynamics
    group : str
        mover group name, default "shooting"
    replace : bool
        whether to replace existing movers in the group; default True
    """
    _level = levels.MOVER
    MoverClass = paths.OneWayShootingMover
    def __init__(self, selector=None, ensembles=None, engine=None,
                 group="shooting", replace=True):
        super(OneWayShootingStrategy, self).__init__(
            ensembles=ensembles, group=group, replace=replace
        )
        if selector is None:
            selector = paths.UniformSelector()
        self.selector = selector
        self.engine = engine

    def make_movers(self, scheme):
        parameters = self.get_parameters(scheme=scheme,
                                         list_parameters=[self.selector],
                                         nonlist_parameters=[self.engine])
        shooters = [
            self.MoverClass(
                ensemble=ens,
                selector=sel,
                engine=eng
            ).named(self.MoverClass.__name__ + " " + ens.name)
            for (ens, sel, eng) in parameters
        ]
        return shooters


class ForwardShootingStrategy(OneWayShootingStrategy):
    """
    Strategy for ForwardShooting only.

    Allows choice of shooting point selector. Useful for e.g. Constrained
    Interface shooting.

    Parameters
    ----------
    selector : :class:`.ShootingPointSelector`
        method used to select shooting point
    ensembles : list of :class:`.Ensemble`
        ensembles for which this strategy applies; None gives default
        behavior
    engine : :class:`.DynamicsEngine`
        engine for the dynamics
    group : str
        mover group name, default "shooting"
    replace : bool
        whether to replace existing movers in the group; default True
    """
    _level = levels.MOVER
    MoverClass = paths.ForwardShootMover


class TwoWayShootingStrategy(SingleEnsembleMoveStrategy):
    """Strategy to make a group of 2-way shooting movers.

    Parameters
    ----------
    modifier : :class:`.SnapshotModifier`
        how to modify the shooting point
    selector : :class:`.ShootingPointSelector`
        how to select the shooting point; None (default) gives uniform
        selection
    ensembles : list of :class:`.Ensemble`
        ensembles to include; see :class:`.MoveStrategy` documentation for
        details
    engine : :class:`.DynamicsEngine`
        the dynamics engine to use
    group : string
        the name of the mover group, default is "shooting"
    replace : bool
        whether to replace existing movers, default True. See
        :class:`.MoveStrategy` documentation for details.
    """
    _level = levels.MOVER
    def __init__(self, modifier, selector=None, ensembles=None, engine=None,
                 group="shooting", replace=True):
        super(TwoWayShootingStrategy, self).__init__(
            ensembles=ensembles, group=group, replace=replace
        )
        self.modifier = modifier
        if selector is None:
            selector = paths.UniformSelector()
        self.selector = selector
        self.engine = engine

    def make_movers(self, scheme):
        parameters = self.get_parameters(
            scheme=scheme,
            list_parameters=[self.selector, self.modifier],
            nonlist_parameters=[self.engine]
        )
        shooters = [
            paths.TwoWayShootingMover(
                ensemble=ens,
                selector=sel,
                modifier=mod,
                engine=eng
            ).named("TwoWayShooting " + ens.name)
            for (ens, sel, mod, eng) in parameters
        ]
        return shooters


class NearestNeighborRepExStrategy(MoveStrategy):
    """
    Make the NN replica exchange scheme among ordered ensembles.
    """
    _level = levels.SIGNATURE
    def __init__(self, ensembles=None, group="repex", replace=True):
        super(NearestNeighborRepExStrategy, self).__init__(
            ensembles=ensembles, group=group, replace=replace
        )

    def make_movers(self, scheme):
        ensemble_list = self.get_ensembles(scheme, self.ensembles)
        movers = []
        for ens in ensemble_list:
            movers.extend([paths.ReplicaExchangeMover(ensemble1=ens[i],
                                                      ensemble2=ens[i+1])
                           for i in range(len(ens)-1)])
        return movers

class NthNearestNeighborRepExStrategy(MoveStrategy):
    _level = levels.SIGNATURE
    pass

# inherits from NearestNeighbor so it can get the same __init__ & _level
class AllSetRepExStrategy(NearestNeighborRepExStrategy):
    """
    Make the replica exchange strategy with all ensembles in each sublist.

    Default is to take a list with each transition (interface set) in a
    different sublist. This makes all the exchanges within that list.
    """
    def make_movers(self, scheme):
        ensemble_list = self.get_ensembles(scheme, self.ensembles)
        movers = []
        for ens in ensemble_list:
            pairs = list(itertools.combinations(ens, 2))
            movers.extend([paths.ReplicaExchangeMover(ensemble1=pair[0],
                                                      ensemble2=pair[1])
                           for pair in pairs])
        return movers

class SelectedPairsRepExStrategy(MoveStrategy):
    """
    Add replica exchange swap for specific pairs of ensembles.

    Note that unlike many signature-level strategies, this defaults to
    `replace=False`, under the assumption that you're probably using it to
    add extra replica exchanges.

    Parameters
    ----------
    ensembles : list of pairs of ensembles, or list of two ensembles
        ensemble pairs for replica exchange
    group : str
        name of the group, default 'repex'
    replace : bool
        whether to replica existing signature, default False
    """
    _level = levels.SIGNATURE

    def initialization_error(self):
        raise RuntimeError("SelectedPairsRepExStrategy must be "
                           + "initialized with ensemble pairs.")

    def __init__(self, ensembles=None, group="repex", replace=False):
        # check that we have a list of pairs
        if ensembles is None:
            self.initialization_error()
        else:
            for pair in ensembles:
                try:
                    pair_len = len(pair)
                except TypeError:
                    pair_len = len(ensembles)
                if pair_len != 2:
                    self.initialization_error()

        super(SelectedPairsRepExStrategy, self).__init__(
            ensembles=ensembles, group=group, replace=replace
        )

    def make_movers(self, scheme):
        ensemble_list = self.get_ensembles(scheme, self.ensembles)
        movers = []
        for pair in ensemble_list:
            movers.append(paths.ReplicaExchangeMover(ensemble1=pair[0],
                                                     ensemble2=pair[1]))
        return movers


class StateSwapRepExStrategy(MoveStrategy):
    pass

class ReplicaExchangeStrategy(MoveStrategy):
    """
    Converts EnsembleHops to ReplicaExchange (single replica to default)
    """
    _level = levels.SUPERGROUP
    def __init__(self, ensembles=None, group="repex", replace=True,
                 from_group=None, bias=None):
        super(ReplicaExchangeStrategy, self).__init__(
            ensembles=ensembles, group=group, replace=replace
        )
        self.bias = bias
        self.from_group = from_group
        if self.from_group is None:
            self.from_group = self.group

    def check_for_hop_repex_validity(self, signatures):
        """
        Checks that the given set of signatures can be either repex or hop.
        """
        # nested function used for error handling
        def sig_error(sig, errstr=""):
            raise RuntimeError("Signature error: " + errstr + str(sig))

        for sig in signatures:
            # We use the fact that Python uses short-circuit logic to throw
            # the exception if the test fails, and the assertion acts as a
            # backup. This is faster than a try: except: version.  see
            # http://stackoverflow.com/questions/1569049/#1569618
            assert(len(sig[0]) == len(sig[1]) or sig_error(sig))
            n_ens = len(sig[0])
            if n_ens == 2:  # replica exchange
                assert(
                    set(sig[0]) == set(sig[1]) or
                    sig_error(sig, errstr="Not replica exchange signature. ")
                )
            elif n_ens == 1: # already ensemble hop (ish)
                assert(
                    # TODO: add test for this
                    (sig[1], sig[0]) in signatures or
                    sig_error(sig, errstr="No detailed balance partner. ")
                )
            else:
                sig_error(sig, errstr="Signature contains " + str(n_ens) + 
                          " ensembles.")

    def make_movers(self, scheme):
        signatures = [m.ensemble_signature
                      for m in scheme.movers[self.from_group]]
        # a KeyError here indicates that there is no existing group of that
        # name: build scheme.movers[self.from_group] before trying to use it!

        self.check_for_hop_repex_validity(signatures)

        swap_list = []
        for sig in signatures:
            n_ens = len(sig[0])
            if n_ens == 2:
                swap_list.append(sig[0])
            elif n_ens == 1:
                swap = (sig[0][0], sig[1][0])
                if not (swap[1], swap[0]) in swap_list:
                    swap_list.append(swap)

        swaps = [paths.ReplicaExchangeMover(swap[0], swap[1])
                 for swap in swap_list]
        return swaps


class EnsembleHopStrategy(ReplicaExchangeStrategy):
    """
    Converts ReplicaExchange to EnsembleHop.

    from_group: can differ from output group `group` if desired
    """
    _level = levels.SUPERGROUP

    def make_movers(self, scheme):
        signatures = [m.ensemble_signature 
                      for m in scheme.movers[self.from_group]]
        # a KeyError here indicates that there is no existing group of that
        # name: build scheme.movers[self.from_group] before trying to use it!

        self.check_for_hop_repex_validity(signatures)

        hop_list = []
        for sig in signatures:
            n_ens = len(sig[0])
            if n_ens == 2:
                hop_list.extend([[sig[0][0],sig[0][1]], [sig[0][1],sig[0][0]]])
            elif n_ens == 1:
                hop_list.extend([[sig[0][0], sig[1][0]]])

        hops = []
        for hop in hop_list:
            if self.bias is not None:
                bias = self.bias.bias_probability(hop[0], hop[1])
            else:
                bias = None
            hopper = paths.EnsembleHopMover(hop[0], hop[1], bias=bias)
            hopper.named("EnsembleHop " + str(hop[0].name) + "->" +
                         str(hop[1].name))
            hops.append(hopper)

        return hops


class PathReversalStrategy(MoveStrategy):
    """
    Creates PathReversalMovers for the strategy.
    """
    _level = levels.MOVER
    def __init__(self, ensembles=None, group="pathreversal", replace=True):
        super(PathReversalStrategy, self).__init__(
            ensembles=ensembles, group=group, replace=replace
        )

    def make_movers(self, scheme):
        ensemble_list = self.get_ensembles(scheme, self.ensembles)
        ensembles = reduce(list.__add__, map(lambda x: list(x), ensemble_list))
        movers = paths.PathReversalSet(ensembles)
        return movers


class MinusMoveStrategy(MoveStrategy):
    """
    Takes a given scheme and makes the minus mover.
    """
    _level = levels.MOVER
    MoverClass = paths.MinusMover
    def __init__(self, engine=None, ensembles=None, group="minus",
                 replace=True):
        super(MinusMoveStrategy, self).__init__(
            ensembles=ensembles, group=group, replace=replace
        )
        self.engine = engine

    def get_ensembles(self, scheme, ensembles):
        network = scheme.network
        if ensembles is None:
            minus_ensembles = network.minus_ensembles
            state_sorted_minus = collections.defaultdict(list)
            for minus in minus_ensembles:
                state_sorted_minus[minus.state_vol].append(minus)
            ensembles = list(state_sorted_minus.values())

        # now we use super's ability to turn it into list-of-list
        res_ensembles = super(MinusMoveStrategy, self).get_ensembles(scheme,
                                                                     ensembles)
        return res_ensembles

    def get_per_mover_ensembles(self, scheme):
        minus_ensembles = self.get_ensembles(scheme, self.ensembles)
        ensembles = reduce(list.__add__,
                           map(lambda x: list(x), minus_ensembles))
        innermosts = []
        special_minus = scheme.network.special_ensembles['minus']
        for ens in ensembles:
            innermosts.append([t.ensembles[0] for t in special_minus[ens]])
        return list(zip(ensembles, innermosts))


    def make_movers(self, scheme):
        parameters = self.get_parameters(scheme,
                                         nonlist_parameters=[self.engine])
        movers = [
            self.MoverClass(minus_ensemble=ens,
                             innermost_ensembles=innermosts,
                             engine=eng)
            for (ens, innermosts, eng) in parameters
        ]
        return movers

class SingleReplicaMinusMoveStrategy(MinusMoveStrategy):
    """
    Takes a given scheme and makes a single-replica minus mover.
    """
    MoverClass = paths.SingleReplicaMinusMover


class OrganizeByMoveGroupStrategy(MoveStrategy):
    """
    Default global strategy.

    First choose move type, then choose specific instance of the mover.

    Attributes
    ----------
    default_group_weights : dict
        In the format {str(group_name) : float(weight)}
    group_weights : dict
        The sortkey weights. In the format {str(group_name) : float(weight)}
    mover_weights = dict
        The mover weights. In the format {(str(group_name),
        ensemble_signature) : weight}
    """
    _level = levels.GLOBAL
    default_group_weights = {
        'shooting' : 1.0,
        'repex' : 0.5,
        'pathreversal' : 0.5,
        'minus' : 0.2
    }
    def __init__(self, ensembles=None, group=None, replace=True):
        super(OrganizeByMoveGroupStrategy, self).__init__(ensembles,
                                                          group, replace)
        self.group_weights = {}
        self.mover_weights = {}
        self.group = group
        self.replace = replace

    def make_chooser(self, scheme, mover_weights, choosername):
        """
        Make RandomChoiceMover based on the movers and weights in
        mover_weights.
        """
        chooser = paths.RandomChoiceMover(
            movers=list(mover_weights.keys()),
            weights=list(mover_weights.values())
        )
        chooser.named(choosername)
        return chooser

    def default_weights(self, scheme):
        """
        Set the default weights given the initial `scheme`.

        Note that this includes preservation of scheme.choice_probability,
        if it is set.

        Parameters
        ----------
        scheme : paths.MoveScheme
            the scheme to which this strategy is being applied

        Returns
        -------
        tuple (sortkey_w, movers_w)
            sortkey_w is a dictionary of sort keys to weights; movers_w is a
            dictionary of mover keys to weights. See class definition for
            the specific formats of the keys.
        """
        sortkey_w = {}
        movers_w = {}
        if scheme.choice_probability != {}:
            # extract weights from the choice probability
            (sortkey_w, movers_w) = self.weights_from_choice_probability(
                scheme, scheme.choice_probability
            )
        else:
            # generate absolutely generic weights
            for skey in scheme.movers.keys():
                sortkey_w[skey] = self.default_group_weights.get(skey, 1.0)
                for mover in scheme.movers[skey]:
                    total_sig = (skey, mover.ensemble_signature)
                    movers_w[total_sig] = 1.0
        return (sortkey_w, movers_w)

    def override_weights(self, weights, override_w):
        """
        Overrides weights in a dictionary.

        TODO: as this got simplified, I think there might be Python
        built-ins to accomplish it (dict.update?)
        """
        for key in override_w.keys():
            weights[key] = override_w[key]
        return weights


    def get_weights(self, scheme, sorted_movers, sort_weights_override=None,
                    mover_weights_override=None):
        """
        Gets sort_weights and mover_weights dictionaries.

        Parameters
        ----------
        scheme : paths.MoveScheme
            the scheme to which this strategy is being applied
        sorted_movers : unneeded?
        sort_weights_override : dict
            Overrides for sort weights. Format {sort_key : weight}; see
            class definition for sort_key format
        mover_weights_override : dict
            Overrides for mover weights. Format {mover_key : weight}; see
            class definition for mover_key format

        Returns
        -------
        tuple (sortkey_w, movers_w)
            Canonical weights for this strategy. sortkey_w is a dictionary
            of sort keys to weights; movers_w is a dictionary of mover keys
            to weights. See class definition for the specific formats of the
            keys.
        """
        if sort_weights_override is None:
            sort_weights_override = dict()

        if mover_weights_override is None:
            mover_weights_override = dict()

        (sorted_w, mover_w) = self.default_weights(scheme)
        sorted_weights = self.override_weights(sorted_w, sort_weights_override)
        mover_weights = self.override_weights(mover_w, mover_weights_override)
        return (sorted_weights, mover_weights) # error if somehow undefined

    def weights_from_choice_probability(self, scheme, choice_probability):
        """Get the contributing weights from existing choice probability.

        Parameters
        ----------
        scheme : paths.MoveScheme
            The scheme to which this strategy is being applied
        choice_probability : dict
            Choice probability dictionary to be separated (typically
            scheme.choice_probability). Format {mover:probability}, where
            probability is normalized over all movers.

        Returns
        -------
        tuple (sortkey_w, movers_w)
            Sort key weights and mover weights consistent with this
            scheme and choice_probability.  sortkey_w is a dictionary of
            sort keys to weights; movers_w is a dictionary of mover keys to
            weights. See class definition for the specific formats of the
            keys.
        """

        # first get the norm-based probabilities, then reset them.
        mover_weights = {}
        group_unscaled = {}
        most_common = {}
        # most_most_common tracks which group has the largest count of
        # common values (used as backup if there is no shooting group)
        for group in scheme.movers:
            group_probs = {m : choice_probability[m]
                           for m in choice_probability
                           if m in scheme.movers[group]}

            # normalize here based on making the most common within the
            # group the baseline (1)
            most_common[group] = most_common_value(list(group_probs.values()))


            for m in group_probs:
                val = group_probs[m] / most_common[group][0]
                mover_weights[(group, m.ensemble_signature)] = val

        for group in scheme.movers:
            m0 = scheme.movers[group][0]
            mover_w0 = mover_weights[(group, m0.ensemble_signature)]
            group_unscaled[group] = choice_probability[m0] / mover_w0

        try:
            scaling = most_common['shooting'][0]
        except KeyError:
            most_most_common = None
            most_most_common_count = 0
            for g in most_common:
                if most_common[g][1] > most_most_common_count:
                    most_most_common_count = most_common[g][1]
                    most_most_common = g
            scaling = most_common[most_most_common][0]

        group_weights = {g : group_unscaled[g] / scaling
                         for g in group_unscaled}

        return (group_weights, mover_weights)

    def choice_probability(self, scheme, group_weights, mover_weights):
        """
        Calculates the probability of choosing to do each move.

        This approach requires that the group_weights and mover_weights
        include all groups and all movers in the actual scheme, otherwise a
        KeyError will occur. This is a safety check. Typically these values
        are obtained from the strategy.get_weights function.
        """
        unnormed = {}
        for g_name in scheme.movers.keys():
            for mover in scheme.movers[g_name]:
                m_sig = (g_name, mover.ensemble_signature)
                unnormed[mover] = group_weights[g_name]*mover_weights[m_sig]
        norm = sum(unnormed.values())
        return {m : unnormed[m] / norm for m in unnormed}

    def chooser_root_weights(self, scheme, group_weights, mover_weights):
        """
        Determine the choice probabilities for the root chooser. The nature
        of the root chooser depends on the class definition.
        """
        weights = {}
        for g in scheme.movers.keys():
            weights[g] = sum([mover_weights[m] for m in mover_weights
                              if m[0] == g]) * group_weights[g]
        return weights

    def chooser_mover_weights(self, scheme, group, mover_weights):
        """
        Set the weights within each "sorted"-level chooser. The nature of
        the sorting depends on the class definition.
        """
        weights = {m : mover_weights[(group, m.ensemble_signature)]
                   for m in scheme.movers[group]}
        return weights

    def make_movers(self, scheme):
        (group_weights, mover_weights) = self.get_weights(
            scheme=scheme,
            sorted_movers=scheme.movers,
            sort_weights_override=self.group_weights,
            mover_weights_override=self.mover_weights
        )
        scheme.choice_probability = self.choice_probability(
            scheme, group_weights, mover_weights
        )
        self.group_weights = group_weights
        self.mover_weights = mover_weights
        root_info = self.chooser_root_weights(scheme, group_weights,
                                              mover_weights)
        chooser_dict = {}
        for group in root_info.keys():
            # care to the order of weights
            weight_dict = self.chooser_mover_weights(scheme, group,
                                                     mover_weights)
            choosername = group.capitalize()+"Chooser"
            chooser_dict[group] = self.make_chooser(scheme, weight_dict,
                                                    choosername)

        root_couples = [(root_info[g], chooser_dict[g])
                        for g in root_info.keys()]
        (root_weights, choosers) = zip(*root_couples)
        root_chooser = paths.RandomChoiceMover(movers=choosers,
                                               weights=root_weights)
        root_chooser.named("RootMover")
        scheme.root_mover = root_chooser
        return root_chooser


class OrganizeByEnsembleStrategy(OrganizeByMoveGroupStrategy):
    """
    Global strategy to organize by ensemble first. Needed for SRTIS.

    First we choose an ensemble, then we choose the specific mover within
    that ensemble.

    Attributes
    ----------
    ensemble_weights : dict
        The sortkey weights. In the format {paths.Ensemble : float(weight)}
    mover_weights : dict
        The mover weights. In the fromat {(str(groupname),
        PathMover.ensemble_signature, paths.Ensemble) : float(weight)}

    """
    def __init__(self, ensembles=None, group=None, replace=True):
        super(OrganizeByEnsembleStrategy, self).__init__(
            ensembles=ensembles, group=group, replace=replace
        )
        self.mover_weights = {}
        self.ensemble_weights = {}

    def weights_from_choice_probability(self, scheme, choice_probability):
        # NOTE this is harder than it looks. The problem is that there isn't
        # always a unique solution when one move appears under more than one
        # ensemble. More details on the problem and solution are in a gist:
        # https://gist.github.com/dwhswenson/5d5b18ba8e811cbe21da
        ensemble_weights = {}
        mover_weights = {}
        ensemble_list = []
        for m in choice_probability:
            ensemble_list.extend([e for e in m.ensemble_signature[0]])
        ensembles = set(ensemble_list)
        for ens in ensembles:
            ens_movers = [m for m in choice_probability
                          if ens in m.ensemble_signature[0]]
            for m in ens_movers:
                ens_sig = m.ensemble_signature
                group = [g for g in scheme.movers if m in scheme.movers[g]][0]
                weight = choice_probability[m] / len(ens_sig[0])
                mover_weights[(group, ens_sig, ens)] = weight

            local_movers = {s : mover_weights[s] for s in mover_weights
                            if s[2] == ens}
            ensemble_weights[ens] = sum(local_movers.values())

            shooters = [s for s in local_movers if s[0] == 'shooting']
            if len(shooters) > 0:
                renorm = local_movers[shooters[0]]
            else:
                renorm = most_common_value(local_movers.values())[0]
            for s in local_movers:
                mover_weights[s] = local_movers[s] / renorm
        ensemble_norm = most_common_value(ensemble_weights.values())[0]
        ensemble_weights = {e : ensemble_weights[e] / ensemble_norm
                            for e in ensemble_weights}
        return (ensemble_weights, mover_weights)

    def default_weights(self, scheme):
        """
        Set the default weights given the initial `scheme`.

        Note that this includes preservation of scheme.choice_probability,
        if it is set.

        Parameters
        ----------
        scheme : paths.MoveScheme
            the scheme to which this strategy is being applied

        Returns
        -------
        tuple (sortkey_w, movers_w)
            sortkey_w is a dictionary of sort keys to weights; movers_w is a
            dictionary of mover keys to weights. See class definition for
            the specific formats of the keys.
        """
        ensemble_weights = {}
        mover_weights = {}
        if scheme.choice_probability != {}:
            (ensemble_weights, mover_weights) = (
                self.weights_from_choice_probability(scheme,
                                                     scheme.choice_probability)
            )
        else:
            for g in scheme.movers:
                for m in scheme.movers[g]:
                    for e in m.ensemble_signature[0]:
                        ensemble_weights[e] = 1.0
                        mover_weights[(g, m.ensemble_signature, e)] = 1.0

        return (ensemble_weights, mover_weights)


    def choice_probability(self, scheme, ensemble_weights, mover_weights):
        """Get the contributing weights from existing choice probability.

        Parameters
        ----------
        scheme : paths.MoveScheme
            The scheme to which this strategy is being applied
        choice_probability : dict
            Choice probability dictionary to be separated (typically
            scheme.choice_probability). Format {mover:probability}, where
            probability is normalized over all movers.

        Returns
        -------
        tuple (sortkey_w, movers_w)
            Sort key weights and mover weights consistent with this
            scheme and choice_probability.  sortkey_w is a dictionary of
            sort keys to weights; movers_w is a dictionary of mover keys to
            weights. See class definition for the specific formats of the
            keys.
        """
        choice_probability = collections.defaultdict(float)
        ens_prob_norm = sum(ensemble_weights.values())
        ens_prob = {e : ensemble_weights[e] / ens_prob_norm
                    for e in ensemble_weights}

        mover_norm = {e : sum([mover_weights[s] for s in mover_weights
                               if s[2] == e])
                      for e in ensemble_weights}

        for sig in mover_weights:
            group, ens_sig, ens = sig
            local_prob = ens_prob[ens] * mover_weights[sig] / mover_norm[ens]
            # take first, because as MacLeod says, "there can be only one!"
            mover = [m for m in scheme.movers[group]
                     if m.ensemble_signature == ens_sig][0]
            choice_probability[mover] += local_prob

        return dict(choice_probability)

    def chooser_root_weights(self, scheme, ensemble_weights, mover_weights):
        """
        Determine the choice probabilities for the root chooser. The nature
        of the root chooser depends on the class definition.
        """
        return ensemble_weights

    def chooser_mover_weights(self, scheme, ensemble, mover_weights):
        """
        Set the weights within each "sorted"-level chooser. The nature of
        the sorting depends on the class definition.
        """
        weights = {}
        for sig in [s for s in mover_weights if s[2] == ensemble]:
            group = sig[0]
            ens_sig = sig[1]
            #ens = sig[2]
            # there can be only one
            mover = [m for m in scheme.movers[group]
                     if m.ensemble_signature == ens_sig][0]
            weights[mover] = mover_weights[sig]
        return weights

    def make_movers(self, scheme):
        (ensemble_weights, mover_weights) = self.get_weights(
            scheme=scheme,
            sorted_movers=scheme.movers,
            sort_weights_override=self.ensemble_weights,
            mover_weights_override=self.mover_weights
        )
        scheme.choice_probability = self.choice_probability(
            scheme, ensemble_weights, mover_weights
        )
        self.ensemble_weights = ensemble_weights
        self.mover_weights = mover_weights
        root_info = self.chooser_root_weights(scheme, ensemble_weights,
                                              mover_weights)

        chooser_dict = {}
        for ens in root_info.keys():
            weight_dict = self.chooser_mover_weights(scheme, ens,
                                                     mover_weights)
            choosername = ens.name+"Chooser"
            chooser_dict[ens] = self.make_chooser(scheme, weight_dict,
                                                  choosername)


        root_couples = [(root_info[g], chooser_dict[g])
                        for g in root_info.keys()]
        (root_weights, choosers) = zip(*root_couples)
        root_chooser = paths.RandomChoiceMover(movers=choosers,
                                               weights=root_weights)
        root_chooser.named("RootMover")
        scheme.root_mover = root_chooser
        return root_chooser

class PoorSingleReplicaStrategy(OrganizeByEnsembleStrategy):
    """
    Organizes by ensemble, then readjusts the weights to have a bunch of
    null moves.
    """
    def __init__(self, ensembles=None, group=None, replace=True):
        super(PoorSingleReplicaStrategy, self).__init__(
            ensembles=ensembles, group=group, replace=replace
        )
        self.null_mover = paths.IdentityPathMover(counts_as_trial=False)

    def chooser_mover_weights(self, scheme, ensemble, mover_weights):
        # this is where I'll have to pad with the null_mover
        if scheme.choice_probability == {}:
            msg = "Set choice probability before chooser_mover_weights"
            raise RuntimeError(msg)

        weights = {}
        for sig in [s for s in mover_weights if s[2] == ensemble]:
            group = sig[0]
            ens_sig = sig[1]
            #ens = sig[2]
            # there can be only one
            mover = [m for m in scheme.movers[group]
                     if m.ensemble_signature == ens_sig][0]
            weights[mover] = scheme.choice_probability[mover]

        sum_weights = sum(weights.values())
        weights[self.null_mover] = 1.0 - sum_weights
        return weights


    def make_movers(self, scheme):
        old_root = super(PoorSingleReplicaStrategy, self).make_movers(scheme)
        root = paths.RandomAllowedChoiceMover(
            movers=old_root.movers,
            weights=old_root.weights
        )
        scheme.real_choice_probability = {
            m : scheme.choice_probability[m] / float(len(root.movers))
            for m in scheme.choice_probability.keys()
        }
        return root


