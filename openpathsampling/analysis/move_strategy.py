import openpathsampling as paths
from openpathsampling.todict import OPSNamed

class MoveStrategy(object):
    def __init__(self, network, scheme=None):
        pass

    def make_chooser(self, scheme, groupname, choosername):
        pass

class UniformSelectionShootingStrategy(MoveStrategy):
    pass

class NearestNeighborRepExStrategy(MoveStrategy):
    pass

class SingleReplicaStrategy(MoveStrategy):
    pass

class ReplicaExchangeStrategy(MoveStrategy):
    pass

class MSOuterRepExStrategy(MoveStrategy):
    pass

class MinusMoveStrategy(MoveStrategy):
    pass

class BiasedSelectionShootingStrategy(MoveStrategy):
    pass

class AllSetRepExStrategy(MoveStrategy):
    pass

