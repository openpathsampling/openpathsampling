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
        return self._level

    @level.setter
    def level(self, value):
        self._level = value
        self.set_replace(self.replace)

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
            input ensembles

        Returns
        -------
        list of list of Ensembles
            regularized output
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
        # TODO: add to hidden ensembles
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
        for weights in [scheme.mover_weights, self.mover_weights]:
            if weights != {}:
                for k in weights.keys():
                    mover_weights[k] = weights[k]

        return mover_weights

    def get_ensemble_weights(self, scheme):
        ensemble_weights = {}
        for group in scheme.movers.keys():
            ensemble_weights[group] = {m.ensemble_signature : 1.0 
                                       for m in scheme.movers[group]} 

        for weights in [scheme.ensemble_weights, self.ensemble_weights]:
            if weights != {}:
                for group in weights.keys():
                    for sig in weights[group].keys():
                        ensemble_weights[group][sig] = weights[group][sig]
        return ensemble_weights



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
            try:
                group_weights = mover_weights[group]
            except KeyError:
                group_weights = 1.0
                mover_weights[group] = group_weights
            weights.append(len(scheme.movers[group])*group_weights)
        root_chooser = paths.RandomChoiceMover(movers=choosers,
                                               weights=weights)
        root_chooser.name = "RootMover"
        scheme.root_mover = root_chooser
        scheme.mover_weights = mover_weights
        scheme.ensemble_weights = ensemble_weights
        return root_chooser


class SingleReplicaStrategy(MoveStrategy):
    """
    Converts ReplicaExchange to EnsembleHop, and changes overall structure
    to SRTIS approach.
    """
    pass

