import openpathsampling as paths
from openpathsampling.todict import OPSNamed
from openpathsampling import PathMoverFactory as pmf

import itertools

import collections
LevelLabels = collections.namedtuple(
    "LevelLabels", 
    ["SIGNATURE", "MOVER", "GROUP", "SUPERGROUP", "GLOBAL"]
)
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
        indices = [i for i in range(len(distances)) if distances[i]==mindist]
        if len(indices) > 1:
            return None
        else:
            return levels[indices[0]]

# possible rename to make available as paths.stategy_levels?
levels = StrategyLevels(
    SIGNATURE=10,
    MOVER=30,
    GROUP=50,
    SUPERGROUP=70,
    GLOBAL=90
)

class MoveStrategy(object):
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
    def __init__(self, ensembles, group, replace, network):
        self.ensembles = ensembles
        self.network = network
        self.group = group
        self.replace = replace
        self.replace_signatures = None
        self.replace_movers = None
        self.set_replace(replace)

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
        self.set_replace(self.replace) # behavior of replace depends on level

    def set_replace(self, replace):
        """Sets values for replace_signatures and replace_movers."""
        self.replace_signatures = False
        self.replace_movers = False
        level_type = levels.level_type(self.level)
        if level_type == levels.SIGNATURE:
            self.replace_signatures = replace
        elif level_type == levels.MOVER:
            self.replace_movers = replace

    def get_ensembles(self, ensembles):
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
            for t in self.network.sampling_transitions:
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
 
class OneWayShootingStrategy(MoveStrategy):
    """
    Strategy for OneWayShooting. Allows choice of shooting point selector.
    """
    _level = levels.MOVER
    def __init__(self, selector=None, ensembles=None, group="shooting",
                 replace=True, network=None):
        super(OneWayShootingStrategy, self).__init__(
            ensembles=ensembles, network=network, group=group, replace=replace
        )
        if selector is None:
            selector = paths.UniformSelector()
        self.selector = selector

    def make_movers(self, scheme):
        if self.network is None:
            self.network = scheme.network
        ensemble_list = self.get_ensembles(self.ensembles)
        ensembles = reduce(list.__add__, map(lambda x: list(x), ensemble_list))
        shooters = pmf.OneWayShootingSet(self.selector, ensembles)
        return shooters

class NearestNeighborRepExStrategy(MoveStrategy):
    """
    Make the NN replica exchange scheme among ordered ensembles.
    """
    _level = levels.SIGNATURE
    def __init__(self, ensembles=None, group="repex", replace=True,
                 network=None):
        super(NearestNeighborRepExStrategy, self).__init__(
            ensembles=ensembles, network=network, group=group, replace=replace
        )

    def make_movers(self, scheme):
        if self.network is None:
            self.network = scheme.network
        ensemble_list = self.get_ensembles(self.ensembles)
        movers = []
        for ens in ensemble_list:
            movers.extend(
                [paths.ReplicaExchangeMover(ensembles=[[ens[i], ens[i+1]]])
                 for i in range(len(ens)-1)]
            )
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
        if self.network is None:
            self.network = scheme.network
        ensemble_list = self.get_ensembles(self.ensembles)
        movers = []
        for ens in ensemble_list:
            pairs = list(itertools.combinations(ens, 2))
            movers.extend([paths.ReplicaExchangeMover(ensembles=list(pair))
                           for pair in pairs])
        return movers

class SelectedPairsRepExStrategy(MoveStrategy):
    """
    Take a specific pair of ensembles and add a replica exchange swap for
    that pair.
    """
    _level = levels.SIGNATURE

    def initialization_error(self):
        raise RuntimeError("SelectedPairsRepExStrategy must be initialized with ensemble pairs.")

    def __init__(self, ensembles=None, group="repex", replace=False,
                 network=None):
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
            ensembles=ensembles, network=network, group=group, replace=replace
        )

    def make_movers(self, scheme):
        if self.network is None:
            self.network = scheme.network
        ensemble_list = self.get_ensembles(self.ensembles)
        movers = []
        for pair in ensemble_list:
            movers.append(paths.ReplicaExchangeMover(ensembles=pair))
        return movers


class StateSwapRepExStrategy(MoveStrategy):
    pass

class ReplicaExchangeStrategy(MoveStrategy):
    _level = levels.SUPERGROUP
    """
    Converts EnsembleHops to ReplicaExchange (single replica to default)
    """
    pass

class EnsembleHopStrategy(MoveStrategy):
    _level = levels.SUPERGROUP
    """
    Converts ReplicaExchange to EnsembleHop.
    """
    pass

class PathReversalStrategy(MoveStrategy):
    """
    Creates PathReversalMovers for the strategy.
    """
    _level = levels.MOVER
    def __init__(self, ensembles=None, group="pathreversal", replace=True,
                 network=None):
        super(PathReversalStrategy, self).__init__(
            ensembles=ensembles, network=network, group=group, replace=replace
        )

    def make_movers(self, scheme):
        if self.network is None:
            self.network = scheme.network
        ensemble_list = self.get_ensembles(self.ensembles)
        ensembles = reduce(list.__add__, map(lambda x: list(x), ensemble_list))
        movers = paths.PathReversalSet(ensembles)
        return movers


class MinusMoveStrategy(MoveStrategy):
    """
    Takes a given scheme and makes the minus mover.
    """
    _level = levels.MOVER
    def __init__(self, ensembles=None, group="minus", replace=True,
                 network=None):
        super(MinusMoveStrategy, self).__init__(
            ensembles=ensembles, network=network, group=group, replace=replace
        )

    def get_ensembles(self, ensembles):
        network = self.network
        if ensembles is None:
            minus_ensembles = network.minus_ensembles
            state_sorted_minus = {}
            for minus in minus_ensembles:
                try:
                    state_sorted_minus[minus.state_vol].append(minus)
                except KeyError:
                    state_sorted_minus[minus.state_vol] = [minus]
            ensembles = state_sorted_minus.values()

        # now we use super's ability to turn it into list-of-list
        res_ensembles = super(MinusMoveStrategy, self).get_ensembles(ensembles)
        return res_ensembles

    def make_movers(self, scheme):
        if self.network is None:
            self.network = scheme.network
        network = self.network
        ensemble_list = self.get_ensembles(self.ensembles)
        ensembles = reduce(list.__add__, map(lambda x: list(x), ensemble_list))
        movers = []
        for ens in ensembles:
            innermosts = [t.ensembles[0] 
                          for t in network.special_ensembles['minus'][ens]]
            movers.append(paths.MinusMover(
                minus_ensemble=ens, 
                innermost_ensembles=innermosts
            ))
        return movers

class SingleReplicaMinusMoveStrategy(MinusMoveStrategy):
    pass

class OrganizeByEnsembleStrategy(MoveStrategy):
    _level = levels.GLOBAL
    def __init__(self, ensembles=None, group=None, replace=True,
                 network=None):
        super(OrganizeByEnsembleStrategy, self).__init__(
            ensembles=ensembles, network=network, group=group, replace=replace
        )
        self.ensemble_weights = None
        self.mover_weights = None

    def make_ensemble_level_chooser(self, scheme, ensemble):
        pass

    def make_movers(self, scheme):
        if self.network is None:
            self.network = scheme.network
        network = self.network
        # TODO
        pass

class DefaultStrategy(MoveStrategy):
    _level = levels.GLOBAL
    default_mover_weights = {
        'shooting' : 1.0,
        'repex' : 0.5,
        'pathreversal' : 0.5,
        'minus' : 0.2,
    }
    def __init__(self, ensembles=None, group=None, replace=True,
                 network=None):
        self.mover_weights = {}
        self.ensemble_weights = {}
        self.group = group
        self.replace = replace
        self.network = network

    def make_chooser(self, scheme, group, weights=None, choosername=None):
        if choosername is None:
            choosername = group.capitalize()+"Chooser"
        chooser = paths.RandomChoiceMover(
            movers=scheme.movers[group],
            weights=weights
        )
        chooser.name = choosername
        return chooser

    def get_mover_weights(self, scheme):
        """
        TODO: BRIEF

        Start with the defaults. Override with scheme.mover_weights, then
        override with self.mover_weights. Return the result without
        modifying the others.

        Parameters
        ----------
        scheme : MoveScheme

        Returns
        -------
        dict
        """
        mover_weights = self.default_mover_weights
        for weights in [self.mover_weights]:
            if weights != {}:
                for k in weights.keys():
                    mover_weights[k] = weights[k]

        return mover_weights

    def get_ensemble_weights(self, scheme):
        ensemble_weights = {}
        for group in scheme.movers.keys():
            ensemble_weights[group] = {m.ensemble_signature : 1.0 
                                       for m in scheme.movers[group]} 

        for weights in [self.ensemble_weights]:
            if weights != {}:
                for group in weights.keys():
                    for sig in weights[group].keys():
                        ensemble_weights[group][sig] = weights[group][sig]
        return ensemble_weights


    def choice_probability(self, scheme, mover_weights, ensemble_weights):
        unnormed = {}
        for groupname in scheme.movers.keys():
            group = scheme.movers[groupname]
            group_w = mover_weights[groupname] 
            sig_weights = ensemble_weights[groupname]
            for mover in group:
                sig_w = sig_weights[mover.ensemble_signature]
                unnormed[mover] = group_w * sig_w
        norm = sum(unnormed.values())
        return {m : unnormed[m] / norm for m in unnormed}

    def make_movers(self, scheme):
        if self.network is None:
            self.network = scheme.network

        mover_weights = self.get_mover_weights(scheme)
        ensemble_weights = self.get_ensemble_weights(scheme)
        choosers = []
        weights = []
        for group in scheme.movers.keys():
            # care to the order of weights
            ens_weights = [ensemble_weights[group][m.ensemble_signature]
                           for m in scheme.movers[group]]
            choosers.append(
                self.make_chooser(scheme, group, weights=ens_weights)
            )
            n_movers = len(scheme.movers[group])
            try:
                group_weights = mover_weights[group]*n_movers
            except KeyError:
                mover_weights[group] = 1.0
                group_weights = mover_weights[group]*n_movers
            weights.append(group_weights)
        root_chooser = paths.RandomChoiceMover(movers=choosers,
                                               weights=weights)
        root_chooser.name = "RootMover"
        scheme.root_mover = root_chooser
        scheme.choice_probability = self.choice_probability(
            scheme, mover_weights, ensemble_weights
        )
        return root_chooser


class SingleReplicaStrategy(MoveStrategy):
    """
    Converts ReplicaExchange to EnsembleHop, and changes overall structure
    to SRTIS approach.
    """
    pass

