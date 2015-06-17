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
    _level = 0
    def __init__(self, ensembles, group, replace, network):
        self.ensembles = ensembles
        self.network = network
        self.group = group
        self.replace = replace
        self.replace_signatures = False
        self.replace_movers = False
        self.set_replace(replace)

    @property
    def level(self):
        return self._level

    @level.setter
    def level(self, value):
        self._level = value
        self.set_replace(self.replace)

    def set_replace(self, replace):
        """Sets values for replace_signatures and replace_movers.  """
        level_type = levels.level_type(self.level)
        if level_type == levels.SIGNATURE:
            self.replace_signatures = replace
        elif level_type == levels.MOVER:
            self.replace_movers = replace

    # TODO: is this really going to be used yet?
    def make_chooser(self, scheme, group, choosername=None):
        if choosername is None:
            choosername = groupname.capitalize()+"Chooser"
        chooser = paths.RandomChoiceMover(movers=scheme.movers[groupname])
        chooser.name = choosername
        scheme.include_movers([chooser], 'choosers', replace=False)

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
                 for i in range(len(ensembles)-1)]
            )
        return movers

    def make_scheme(self, scheme=None, ensembles=None):
        """
        Make the NN replica exchange scheme among ordered ensembles.

        Parameters
        ----------
        scheme : MoveScheme (None)
            scheme to start from. If `None`, use self.scheme
        ensembles : list of Ensemble (None)
            ordered list of the ensembles; replica exchange moves are made
            for each pair of neighbors in the list. In None, defaults to
            using per-transition ensembles sets from `self.network`
        group : string ("repex")
            name of the mover group for this. 
        replace : bool (False)
            Whether to replace the existing mover group of this name. If
            False, appends moves to the existing group.
        chooser : bool (True)
            Whether to create a default chooser (RandomChoiceMover) for this
            group. The name of the chooser would be "GroupnameChooser".

        Returns
        -------
        MoveScheme :
            the resulting MoveScheme
        """
        if scheme is not None and self.network is None:
            self.network = scheme.network
        movers = self.make_movers()
        scheme.include_movers(movers, groupname, replace)
        if chooser:
            make_chooser(scheme, groupname, choosername)
        return scheme

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
    def __init__(self, ensembles=None, network=None):
        shooting = OneWayShootingStrategy(
            network=network
        )
        pass



class SingleReplicaStrategy(MoveStrategy):
    """
    Converts ReplicaExchange to EnsembleHop, and changes overall structure
    to SRTIS approach.
    """
    pass

