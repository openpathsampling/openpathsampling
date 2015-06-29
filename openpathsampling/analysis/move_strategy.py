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


class DefaultStrategy(MoveStrategy):
    _level = levels.GLOBAL
    default_group_weights = {
        'shooting' : 1.0,
        'repex' : 0.5,
        'pathreversal' : 0.5,
        'minus' : 0.2
    }
    def __init__(self, ensembles=None, group=None, replace=True,
                 network=None):
        self.group_weights = {}
        self.mover_weights = {}
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


    def strategy_group_weights(self, scheme, mover_weights=None):
        """
        Returns the group weights given by info in the strategy and scheme.

        `mover_weights` is None if either (a) scheme.choice_probability is
        not set; or (b) self.group_weights is set. Then we don't make use of
        scheme.choice_probability. If given mover_weights, we use
        scheme.choice_probability to set the group_weights.

        See also
        --------
        DefaultStrategy.get_weights
        DefaultStrategy.strategy_mover_weights
        """
        group_weights = {}
        if mover_weights is None:
            for group in scheme.movers:
                try:
                    group_weights[group] = self.default_group_weights[group]
                except KeyError:
                    group_weights[group] = 1.0
                # override default
                if group in self.group_weights:
                    group_weights[group] = self.group_weights[group]
        else:
            mover_weights = self.strategy_mover_weights(scheme)
            group_weights = {}
            for group in scheme.movers:
                movers = scheme.movers[group]
                group_weights[group] = sum([scheme.choice_probability[m] 
                                            for m in movers])
            try:
                normalizer = group_weights['shooting']
            except KeyError:
                normalizer = sum(group_weights.values())
            for group in group_weights:
                group_weights[group] /= normalizer

        return group_weights

    def strategy_mover_weights(self, scheme, group_weights=None):
        """
        Returns the mover weights given by info in the strategy and scheme.

        `group_weights` is None if either (a) scheme.choice_probability is
        not set; or (b) self.mover_weights is set. Then we don't make use of
        scheme.choice_probability. If given group_weights, we use
        scheme.choice_probability to set the mover_weights.

        See also
        --------
        DefaultStrategy.get_weights
        DefaultStrategy.strategy_group_weights
        """
        mover_weights = {}
        if group_weights is None:
            for group in scheme.movers:
                mover_weights[group] = {}
                movers = scheme.movers[group]
                # set default
                mover_weights[group] = {m.ensemble_signature : 1.0
                                        for m in scheme.movers[group]}
                if group in self.mover_weights:
                    for sig in self.mover_weights[group]:
                        mover_weights[group][sig] = self.mover_weights[group][sig]
        else:
            m_weights = self.strategy_mover_weights(scheme) # defaults
            pred_choice = self.choice_probability(scheme, group_weights,
                                                  m_weights)
            scheme_choice = scheme.choice_probability
            for group in scheme.movers:
                mover_weights[group] = {
                    m.ensemble_signature : scheme_choice[m] / pred_choice[m]
                    for m in scheme.movers[group]
                }

        for group in mover_weights:
            group_min = min(mover_weights[group].values())
            mover_weights[group] = {s : mover_weights[group][s] / group_min
                                    for s in mover_weights[group]}
                                        
        return mover_weights

    def get_weights(self, scheme):
        """
        BRIEF

        Notes
        -----
        The group_weight gives the relative probability of choosing a group;
        the mover_weight gives the relative probability of choosing a mover
        *within* its group.

        Note that only the variables returned from this tell the full story.
        The defaults and the self.* version may not contain the real set of
        groups/movers.
        """
        choice_prob_set = (scheme.choice_probability != {})
        group_set = (self.group_weights != {})
        mover_set = (self.mover_weights != {})
        if (group_set and mover_set) or not choice_prob_set:
            group_weights = self.strategy_group_weights(scheme)
            mover_weights = self.strategy_mover_weights(scheme)
        elif group_set: #choice_prob is set; mover is not set
            # use group_weights & choice_probability to set mover_weights
            group_weights = self.strategy_group_weights(scheme)
            mover_weights = self.strategy_mover_weights(scheme, group_weights)
        elif mover_set: #choice_prob is set; group is not set
            # use mover_weights & choice_probability to set group_weights
            mover_weights = self.strategy_mover_weights(scheme)
            group_weights = self.strategy_group_weights(scheme, mover_weights)
        else: #choice_prob is set, neither group nor mover is set
            # use default mover weights to get the group weights, then use
            # that to get the actual correct mover weights
            m_weights = self.strategy_mover_weights(scheme)
            group_weights = self.strategy_group_weights(scheme, m_weights)
            mover_weights = self.strategy_mover_weights(scheme, group_weights)

        return (group_weights, mover_weights) # error if somehow undefined


    def choice_probability(self, scheme, group_weights, mover_weights):
        unnormed = {}
        group_norm = sum(group_weights.values())
        for groupname in scheme.movers:
            group_w = group_weights[groupname] / group_norm
            sig_weights = mover_weights[groupname]
            sig_norm = sum(mover_weights[groupname].values())
            for mover in scheme.movers[groupname]:
                sig_w = sig_weights[mover.ensemble_signature] / sig_norm
                unnormed[mover] = group_w * sig_w
        norm = sum(unnormed.values())
        return {m : unnormed[m] / norm for m in unnormed}

    def make_movers(self, scheme):
        if self.network is None:
            self.network = scheme.network

        (group_weights, mover_weights) = self.get_weights(scheme)
        choosers = []
        for group in scheme.movers.keys():
            # care to the order of weights
            ens_weights = [mover_weights[group][m.ensemble_signature]
                           for m in scheme.movers[group]]
            choosers.append(
                self.make_chooser(scheme, group, weights=ens_weights)
            )

        root_weights = [group_weights[group] * len(scheme.movers[group])
                        for group in scheme.movers]
        root_chooser = paths.RandomChoiceMover(movers=choosers,
                                               weights=root_weights)
        root_chooser.name = "RootMover"
        scheme.root_mover = root_chooser
        scheme.choice_probability = self.choice_probability(
            scheme, group_weights, mover_weights
        )
        return root_chooser


class SingleReplicaStrategy(MoveStrategy):
    """
    Converts ReplicaExchange to EnsembleHop, and changes overall structure
    to SRTIS approach.
    """
    pass

class OrganizeByEnsembleStrategy(DefaultStrategy):
    def __init__(self, ensembles=None, group=None, replace=True,
                 network=None):
        super(OrganizeByEnsembleStrategy, self).__init__(
            ensembles=ensembles, network=network, group=group, replace=replace
        )
        self.mover_weights = {}
        self.group_weights = {}

    def make_chooser(self, scheme, ensemble, weights=None, choosername=None):
        if choosername is None:
            choosername = ensemble.name.capitalize() + "Chooser"



    def make_movers(self, scheme):
        if self.network is None:
            self.network = scheme.network
        network = self.network
        (group_weights, mover_weights) = self.get_weights(scheme)
        choosers = []
        for ens in scheme.network.all_ensembles():
            pass

        # TODO
        pass
