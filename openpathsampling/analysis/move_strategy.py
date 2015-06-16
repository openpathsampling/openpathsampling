import openpathsampling as paths
from openpathsampling.todict import OPSNamed

MOVERLEVEL = 10
GROUPLEVEL = 50
SUPERGROUPLEVEL = 75
GLOBALLEVEL = 100

class MoveStrategy(object):
    level = "undefined"
    def __init__(self, group, replace, network):
        self.network = network
        self.group = group
        self.replace = replace

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
    level = MOVERLEVEL
    def __init__(self, selector=None, ensembles=None, group="shooting", replace=True, network=None):
        super(OneWayShootingStrategy, self).__init__(
            network=network, group=group, replace=replace
        )
        if selector is None:
            selector = paths.UniformSelector()
        self.selector = selector
        self.ensembles = ensembles

    def make_scheme(self, scheme=None):
        ensemble_list = self.get_ensembles(self.ensembles)
        ensembles = reduce(list.__add__, map(lambda x: list(x), ensemble_list))
        shooters = paths.PathMoverFactory.OneWayShootingSet(self.selector, ensembles)
        scheme.include_movers(shooters, group, replace)
        return scheme

class NearestNeighborRepExStrategy(MoveStrategy):
    level = GROUPLEVEL
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
    level=GROUPLEVEL
    pass

class AllSetRepExStrategy(MoveStrategy):
    level=GROUPLEVEL
    pass

class SelectedPairsRepExStrategy(MoveStrategy):
    level=GROUPLEVEL
    pass

class StateSwapRepExStrategy(MoveStrategy):
    pass

class ReplicaExchangeStrategy(MoveStrategy):
    level=SUPERGROUPLEVEL
    """
    Converts EnsembleHops to ReplicaExchange (single replica to default)
    """
    pass

class EnsembleHopStrategy(MoveStrategy):
    level=SUPERGROUPLEVEL
    """
    Converts ReplicaExchange to EnsembleHop.
    """
    pass

class PathReversalStrategy(MoveStrategy):
    level=GROUPLEVEL
    pass


class MinusMoveStrategy(MoveStrategy):
    """
    Takes a given network and makes the minus mover.
    """
    pass

class SingleReplicaMinusMoveStrategy(MoveStrategy):
    pass

class DefaultStrategy(MoveStrategy):
    level = GLOBALLEVEL
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

