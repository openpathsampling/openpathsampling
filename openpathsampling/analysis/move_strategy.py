import openpathsampling as paths
from openpathsampling.todict import OPSNamed
from openpathsampling import PathMoverFactory as pmf

import collections
LevelLabels = collections.namedtuple(
    "LevelLabels", 
    ["MOVER", "MOVER_GROUP_EDGE", "GROUP", "GROUP_SUPERGROUP_EDGE", 
     "SUPERGROUP", "SUPERGROUP_GLOBAL_EDGE", "GLOBAL"]
)
levels = LevelLabels(
    MOVER=10,
    MOVER_GROUP_EDGE=40,
    GROUP=50,
    GROUP_SUPERGROUP_EDGE=60,
    SUPERGROUP=75,
    SUPERGROUP_GLOBAL_EDGE=90,
    GLOBAL=100
)

class MoveStrategy(object):
    level = "undefined"
    def __init__(self, group, replace, network):
        self.network = network
        self.group = group
        self.replace = replace

    def set_replace(self, replace):
        if self.level < levels.MOVER_GROUP_EDGE:
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
    level = levels.MOVER
    def __init__(self, selector=None, ensembles=None, group="shooting", replace=True, network=None):
        super(OneWayShootingStrategy, self).__init__(
            network=network, group=group, replace=replace
        )
        if selector is None:
            selector = paths.UniformSelector()
        self.selector = selector
        self.ensembles = ensembles

    def make_movers(self, ensembles):
        ensembles = reduce(list.__add__, map(lambda x: list(x), ensembles))
        shooters = pmf.OneWayShootingSet(self.selector, ensembles)
        return shooters


    def make_scheme(self, scheme=None):
        ensemble_list = self.get_ensembles(self.ensembles)
        shooters = self.make_movers(ensembles)
        scheme.include_movers(shooters, self.group, self.replace)
        return scheme

class NearestNeighborRepExStrategy(MoveStrategy):
    level = levels.GROUP
    def __init__(self, group="repex", replace=True, network=None):
        super(NearestNeighborRepExStrategy, self).__init__(
            network=network, group=group, replace=replace
        )

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
        ensemble_list = self.get_ensembles(ensembles)
        movers = []
        for ens in ensemble_list:
            movers.extend(
                [paths.ReplicaExchangeMover(ensembles=[[ens[i], ens[i+1]]])
                 for i in range(len(ensembles)-1)]
            )

        scheme.include_movers(movers, groupname, replace)
        if chooser:
            make_chooser(scheme, groupname, choosername)
        return scheme

class NthNearestNeighborRepExStrategy(MoveStrategy):
    level=levels.GROUP
    pass

class AllSetRepExStrategy(MoveStrategy):
    level=levels.GROUP
    pass

class SelectedPairsRepExStrategy(MoveStrategy):
    level=levels.GROUP
    pass

class StateSwapRepExStrategy(MoveStrategy):
    pass

class ReplicaExchangeStrategy(MoveStrategy):
    level=levels.SUPERGROUP
    """
    Converts EnsembleHops to ReplicaExchange (single replica to default)
    """
    pass

class EnsembleHopStrategy(MoveStrategy):
    level=levels.SUPERGROUP
    """
    Converts ReplicaExchange to EnsembleHop.
    """
    pass

class PathReversalStrategy(MoveStrategy):
    level=levels.GROUP
    pass


class MinusMoveStrategy(MoveStrategy):
    """
    Takes a given network and makes the minus mover.
    """
    pass

class SingleReplicaMinusMoveStrategy(MoveStrategy):
    pass

class DefaultStrategy(MoveStrategy):
    level = levels.GLOBAL
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

