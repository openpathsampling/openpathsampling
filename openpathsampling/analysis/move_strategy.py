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

# possible rename to make available as paths.strategy_levels?
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

    def normalization_basis(self, scheme):
        return scheme.movers['shooting'][0]


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
            print "pre-norm: ", sorted_weights
            sorted_weights = self.normalize_sorted_weights(scheme, 
                                                           sorted_movers,
                                                           sorted_weights)

        return sorted_weights

    def normalize_sorted_weights(self, scheme, sorted_movers, sorted_weights):
        try:
            normalization_basis = self.normalization_basis(scheme)
        except KeyError:
            normalizer = sum(sorted_weights.values())
        else:
            for skey in sorted_movers:
                if normalization_basis in sorted_movers[skey]:
                    normalizer = sorted_weights[skey]

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

        mover_weights = self.normalize_mover_weights(scheme, sorted_movers,
                                                     mover_weights)

        return mover_weights


    def normalize_mover_weights(self, scheme, sorted_movers, mover_weights):
        for skey in mover_weights:
            default_norm = mover_weights[skey].values()[0]
            try:
                normalization_basis = self.normalization_basis(scheme)
            except KeyError:
                normalization = default_norm
            else:
                norm_key = self._mover_key(normalization_basis, scheme)
                if normalization_basis in sorted_movers[skey]:
                    normalization = mover_weights[skey][norm_key]
                else:
                    normalization = default_norm

            mover_weights[skey] = {
                s : mover_weights[skey][s] / normalization
                for s in mover_weights[skey]
            }
        return mover_weights

    def default_weights(self, scheme):
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
                try:
                    sortkey_w[skey] = self.default_group_weights[skey]
                except KeyError:
                    sortkey_w[skey] = 1.0
                for mover in scheme.movers[skey]:
                    total_sig = (skey, mover.ensemble_signature)
                    movers_w[total_sig] = 1.0
        return (sortkey_w, movers_w)

    def override_weights(self, weights, override_w):
        for key in override_w.keys():
            weights[key] = override_w[key]
        return weights


    def get_weights(self, scheme, sorted_movers, sort_weights_override={}, 
                    mover_weights_override={}):
        """
        Gets sort_weights and mover_weights dictionaries.

        Notes
        -----
        Note that only the variables returned from this tell the full story.
        The defaults and the self.* version may not contain the real set of
        groups/movers.
        """
        (sorted_w, mover_w) = self.default_weights(scheme)
        sorted_weights = self.override_weights(sorted_w, sort_weights_override)
        mover_weights = self.override_weights(mover_w, mover_weights_override)
        return (sorted_weights, mover_weights) # error if somehow undefined

    def weights_from_choice_probability(self, scheme, choice_probability):
        """Get the contributing weights from existing choice probability
        """
        # first get the norm-based probabilities, then reset them.
        mover_weights = {}
        group_unscaled = {}
        most_common = {}
        # most_most_common tracks which group has the largest count of
        # common values (used as backup if there is no shooting group)
        most_most_common = None
        most_most_common_count = None
        for group in scheme.movers:
            group_probs = {m : choice_probability[m] 
                           for m in choice_probability 
                           if m in scheme.movers[group]}
            # normalize here based on making the most common within the
            # group the baseline (1)
            counts = {}
            for v in group_probs.values():
                try:
                    counts[v] += 1
                except:
                    counts[v] = 1
            most_common[group] = None
            most_common_count = 0
            for v in counts:
                if counts[v] > most_common_count:
                    most_common[group] = v
                    if counts[v] > most_most_common_count:
                        most_most_common = group
                        most_most_common_count = counts[v]

            for m in group_probs:
                val = group_probs[m] / most_common[group]
                mover_weights[(group, m.ensemble_signature)] = val

        for group in scheme.movers:
            m0 = scheme.movers[group][0]
            mover_w0 = mover_weights[(group, m0.ensemble_signature)]
            group_unscaled[group] = choice_probability[m0] / mover_w0 

        try:
            scaling = most_common['shooting']
        except KeyError:
            scaling = most_common[most_most_common]

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
        Determine the choice probabilities for the root chooser.

        In this case, the root chooser selects move type.
        """
        weights = {}
        for g in scheme.movers.keys():
            weights[g] = sum([mover_weights[m] for m in mover_weights 
                              if m[0]==g])  * group_weights[g]
        return weights

    def chooser_mover_weights(self, scheme, group, mover_weights):
        weights = {m : mover_weights[(group, m.ensemble_signature)]
                   for m in scheme.movers[group]}
        return weights

    def make_movers(self, scheme):
        (group_weights, mover_weights) = self.get_weights(
            scheme=scheme, 
            sorted_movers=scheme.movers, 
            sort_weights_override=self.group_weights
        )
        scheme.choice_probability = self.choice_probability(
            scheme, group_weights, mover_weights
        )
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
        root_chooser.name = "RootMover"
        scheme.root_mover = root_chooser
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

    # TODO: make this the same for both (this is unique) and maybe move it
    # to scheme.mover_key
    def _mover_key(self, mover, scheme):
        # take first, because as MacLeod says, "there can be only one!"
        group = [g for g in scheme.movers if mover in scheme.movers[g]][0]
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
            sorted_movers=ensemble_movers, 
            sort_weights_override=self.ensemble_weights
        )
        ens_list = mover_weights.keys() # used for canonical ordering
        # (otherwise the weight list isn't in the same order as choosers!)
        for ens in ens_list:
            for mover_key in mover_weights[ens]:
                n_inp = len(mover_key[1][0])
                mover_weights[ens][mover_key] /= n_inp

        choosers = []
        for ens in ens_list:
            weight_dict = {m : mover_weights[ens][self._mover_key(m, scheme)]
                           for m in ensemble_movers[ens]}
            choosername = ens.name.capitalize() + " Chooser"
            choosers.append(
                self.make_chooser(scheme, weight_dict, choosername)
            )
        
        # this sum, among other things, ensure that the prob of selecting an
        # ensemble hop stays the same
        corrected_ensemble_weights = {
            e : ensemble_weights[e] * sum(mover_weights[e].values())
            for e in ens_list
        }
        root_weights = [corrected_ensemble_weights[e] for e in ens_list]
        root_chooser = paths.RandomChoiceMover(movers=choosers,
                                               weights=root_weights)
        root_chooser.name = "RootMover"
        scheme.root_mover = root_chooser
        scheme.choice_probability = self.choice_probability(
            scheme, ensemble_movers, corrected_ensemble_weights, mover_weights
        )
        return root_chooser

