import openpathsampling as paths
from openpathsampling.todict import OPSNamed
from openpathsampling import PathMoverFactory as pmf

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
            for elem in ensembles:
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

class AllSetRepExStrategy(MoveStrategy):
    _level = levels.SIGNATURE
    pass

class SelectedPairsRepExStrategy(MoveStrategy):
    _level = levels.SIGNATURE
    pass

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
    _level = levels.SIGNATURE
    pass


class MinusMoveStrategy(MoveStrategy):
    """
    Takes a given network and makes the minus mover.
    """
    pass

class SingleReplicaMinusMoveStrategy(MoveStrategy):
    pass

class DefaultStrategy(MoveStrategy):
    _level = levels.GLOBAL
    def __init__(self, ensembles=None, group=None, replace=True,
                 network=None):
        self.weight_adjustment = {
            'shooting' : 1.0,
            'repex' : 0.5,
            'pathreversal' : 0.5,
            'minus' : 0.2,
        }
        self.group = group
        self.replace = replace
        self.network = network

    def make_chooser(self, scheme, group, choosername=None):
        if choosername is None:
            choosername = group.capitalize()+"Chooser"
        chooser = paths.RandomChoiceMover(movers=scheme.movers[group])
        chooser.name = choosername
        return chooser

    def make_movers(self, scheme):
        if self.network is None:
            self.network = scheme.network
        choosers = []
        weights = []
        for group in scheme.movers.keys():
            choosers.append(self.make_chooser(scheme, group))
            try:
                weight_adjustment = self.weight_adjustment[group]
            except KeyError:
                weight_adjustment = 1.0
            weights.append(len(scheme.movers[group])*weight_adjustment)
        root_chooser = paths.RandomChoiceMover(movers=choosers,
                                               weights=weights)
        root_chooser.name = "RootMover"
        scheme.root_mover = root_chooser
        return root_chooser


class SingleReplicaStrategy(MoveStrategy):
    """
    Converts ReplicaExchange to EnsembleHop, and changes overall structure
    to SRTIS approach.
    """
    pass

