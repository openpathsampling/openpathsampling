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
    def __init__(self, ensembles, group, replace):
        self.ensembles = ensembles
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
 
class OneWayShootingStrategy(MoveStrategy):
    """
    Strategy for OneWayShooting. Allows choice of shooting point selector.
    """
    _level = levels.MOVER
    def __init__(self, selector=None, ensembles=None, group="shooting",
                 replace=True):
        super(OneWayShootingStrategy, self).__init__(
            ensembles=ensembles, group=group, replace=replace
        )
        if selector is None:
            selector = paths.UniformSelector()
        self.selector = selector

    def make_movers(self, scheme):
        ensemble_list = self.get_ensembles(scheme, self.ensembles)
        ensembles = reduce(list.__add__, map(lambda x: list(x), ensemble_list))
        shooters = pmf.OneWayShootingSet(self.selector, ensembles)
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
        ensemble_list = self.get_ensembles(scheme, self.ensembles)
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
    def __init__(self, ensembles=None, group="minus", replace=True):
        super(MinusMoveStrategy, self).__init__(
            ensembles=ensembles, group=group, replace=replace
        )

    def get_ensembles(self, scheme, ensembles):
        network = scheme.network
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
        res_ensembles = super(MinusMoveStrategy, self).get_ensembles(scheme,
                                                                     ensembles)
        return res_ensembles

    def make_movers(self, scheme):
        network = scheme.network
        ensemble_list = self.get_ensembles(scheme, self.ensembles)
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
    """
    Default global strategy. 
    
    First choose move type, then choose specific instance of the mover.
    """
    _level = levels.GLOBAL
    default_group_weights = {
        'shooting' : 1.0,
        'repex' : 0.5,
        'pathreversal' : 0.5,
        'minus' : 0.2
    }
    def __init__(self, ensembles=None, group=None, replace=True):
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
            movers=mover_weights.keys(),
            weights=mover_weights.values()
        )
        chooser.name = choosername
        return chooser


    def strategy_sortkey_weights(self, scheme, sorted_movers,
                                 sortkey_weights, mover_weights=None):
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
        sorted_weights = {}
        if mover_weights is None:
            for skey in sorted_movers:
                try:
                    # note that since ONLY the keys in group_weight have
                    # defaults that are other than 1, this works REGARDLESS
                    # of the choice of sortkey!
                    sorted_weights[skey] = self.default_group_weights[skey]
                except KeyError:
                    # note that this also works if the key is not a group
                    # (e.g., if you organize by ensemble!)
                    sorted_weights[skey] = 1.0
                # override default
                if skey in sortkey_weights:
                    sorted_weights[skey] = sortkey_weights[skey]
        else:
            mover_weights = self.strategy_mover_weights(scheme, sorted_movers)
            sorted_weights = {}
            for skey in sorted_movers:
                movers = sorted_movers[skey]
                sorted_weights[skey] = sum([scheme.choice_probability[m] 
                                            for m in movers])
            try:
                normalizer = sorted_weights['shooting']
            except KeyError:
                normalizer = sum(sorted_weights.values())
            for skey in sorted_weights:
                sorted_weights[skey] /= normalizer

        return sorted_weights

    def _mover_key(self, mover, scheme):
        return mover.ensemble_signature

    def strategy_mover_weights(self, scheme, sorted_movers, 
                               sortkey_weights=None):
        """
        Returns the mover weights given by info in the strategy and scheme.

        `sortkey_weights` is None if either (a) scheme.choice_probability is
        not set; or (b) self.mover_weights is set. Then we don't make use of
        scheme.choice_probability. If given sortkey_weights, we use
        scheme.choice_probability to set the mover_weights.

        See also
        --------
        DefaultStrategy.get_weights
        DefaultStrategy.strategy_sortkey_weights
        """
        mover_weights = {}
        if sortkey_weights is None:
            for sortkey in sorted_movers:
                mover_weights[sortkey] = {}
                movers = sorted_movers[sortkey]
                # set default
                mover_weights[sortkey] = {self._mover_key(m, scheme): 1.0
                                          for m in sorted_movers[sortkey]}
                if sortkey in self.mover_weights:
                    for sig in self.mover_weights[sortkey]:
                        mover_weights[sortkey][sig] = self.mover_weights[sortkey][sig]
        else:
            m_weights = self.strategy_mover_weights(scheme, sorted_movers) # defaults
            pred_choice = self.choice_probability(scheme, sorted_movers,
                                                  sortkey_weights, m_weights)
            scheme_choice = scheme.choice_probability
            for sortkey in sorted_movers:
                mover_weights[sortkey] = {
                    self._mover_key(m, scheme) : scheme_choice[m]/pred_choice[m]
                    for m in sorted_movers[sortkey]
                }

        for skey in mover_weights:
            skey_min = min(mover_weights[skey].values())
            mover_weights[skey] = {s : mover_weights[skey][s] / skey_min
                                   for s in mover_weights[skey]}
                                        
        return mover_weights

    def get_weights(self, scheme, sorted_movers, preset_sortkey_weights):
        """
        Gets sort_weights and mover_weights dictionaries.

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
        sort_set = (preset_sortkey_weights != {})
        mover_set = (self.mover_weights != {})
        if (sort_set and mover_set) or not choice_prob_set:
            sorted_weights = self.strategy_sortkey_weights(
                scheme, sorted_movers, preset_sortkey_weights
            )
            mover_weights = self.strategy_mover_weights(scheme,
                                                        sorted_movers)
        elif sort_set: #choice_prob is set; mover is not set
            # use sorted_weights & choice_probability to set mover_weights
            sorted_weights = self.strategy_sortkey_weights(
                scheme, sorted_movers, preset_sortkey_weights
            )
            mover_weights = self.strategy_mover_weights(scheme,
                                                        sorted_movers, 
                                                        sorted_weights)
        elif mover_set: #choice_prob is set; sort is not set
            # use mover_weights & choice_probability to set group_weights
            mover_weights = self.strategy_mover_weights(scheme,
                                                        sorted_movers)
            sorted_weights = self.strategy_sortkey_weights(
                scheme, sorted_movers, preset_sortkey_weights, mover_weights
            )
        else: #choice_prob is set, neither group nor mover is set
            # use default mover weights to get the group weights, then use
            # that to get the actual correct mover weights
            m_weights = self.strategy_mover_weights(scheme, sorted_movers)
            sorted_weights = self.strategy_sortkey_weights(
                scheme, sorted_movers, preset_sortkey_weights, m_weights
            )
            mover_weights = self.strategy_mover_weights(scheme,
                                                        sorted_movers, 
                                                        sorted_weights)

        return (sorted_weights, mover_weights) # error if somehow undefined


    def choice_probability(self, scheme, sorted_movers, sorted_weights, 
                           mover_weights):
        """
        Calculates the probability of choosing to do each move.
        """
        unnormed = {}
        sorted_norm = sum(sorted_weights.values())
        for sortkey in sorted_movers:
            sorted_w = sorted_weights[sortkey] / sorted_norm
            sig_weights = mover_weights[sortkey]
            sig_norm = sum(mover_weights[sortkey].values())
            for mover in sorted_movers[sortkey]:
                sig_w = sig_weights[self._mover_key(mover, scheme)] / sig_norm
                unnormed[mover] = sorted_w * sig_w
        norm = sum(unnormed.values())
        return {m : unnormed[m] / norm for m in unnormed}

    def make_movers(self, scheme):
        (group_weights, mover_weights) = self.get_weights(
            scheme=scheme, 
            sorted_movers=scheme.movers, 
            preset_sortkey_weights=self.group_weights
        )
        choosers = []
        for group in scheme.movers.keys():
            # care to the order of weights
            ens_weights = [mover_weights[group][self._mover_key(m, scheme)]
                           for m in scheme.movers[group]]
            weight_dict = {m : w for (m, w) in zip(scheme.movers[group],
                                                   ens_weights)}
            choosername = group.capitalize()+"Chooser"
            choosers.append(
                self.make_chooser(scheme, weight_dict, choosername)
            )

        root_weights = [group_weights[group] * len(scheme.movers[group])
                        for group in scheme.movers]
        root_chooser = paths.RandomChoiceMover(movers=choosers,
                                               weights=root_weights)
        root_chooser.name = "RootMover"
        scheme.root_mover = root_chooser
        scheme.choice_probability = self.choice_probability(
            scheme, scheme.movers, group_weights, mover_weights
        )
        return root_chooser


class SingleReplicaStrategy(MoveStrategy):
    """
    Converts ReplicaExchange to EnsembleHop, and changes overall structure
    to SRTIS approach.
    """
    pass

class OrganizeByEnsembleStrategy(DefaultStrategy):
    def __init__(self, ensembles=None, group=None, replace=True):
        super(OrganizeByEnsembleStrategy, self).__init__(
            ensembles=ensembles, group=group, replace=replace
        )
        self.mover_weights = {}
        self.ensemble_weights = {}

    def _mover_key(self, mover, scheme):
        # take first, because as MacLeod says, "there can be only one!"
        try:
            group = [g for g in scheme.movers if mover in scheme.movers[g]][0]
        except IndexError:
            print mover

        return (group, mover.ensemble_signature)

    def make_movers(self, scheme):
        all_movers = []
        for g in scheme.movers:
            all_movers.extend(scheme.movers[g])

        # ensemble_movers is used in the same way as scheme.movers in the
        # standard organization strategy
        ensemble_movers = {}
        for mover in all_movers:
            for inp_ens in mover.input_ensembles:
                try:
                    ensemble_movers[inp_ens].append(mover)
                except KeyError:
                    ensemble_movers[inp_ens] = [mover]

        (ensemble_weights, mover_weights) = self.get_weights(
            scheme=scheme,
            sorted_movers=ensemble_movers, # TODO: this will change
            preset_sortkey_weights=self.ensemble_weights
        )
        for ens in mover_weights:
            for mover_key in mover_weights[ens]:
                n_inp = len(mover_key[1][0])
                mover_weights[ens][mover_key] /= n_inp
        for ens in mover_weights:
            print "====", ens.name
            print mover_weights[ens]
        #return root
