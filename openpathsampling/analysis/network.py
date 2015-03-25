
class TransitionNetwork(object):
    def __init__(self):
        pass

    def replica_flow_map(self):
        pass

    def replica_flow(self, bottom, top):
        pass


class TISNetwork(TransitionNetwork):
    def __init__(self):
        pass

    def from_transitions(self, transitions, interfaces=None):
        pass

    def default_mover(self):
        pass


class MSTISNetwork(TISNetwork):
    def __init__(self):
        pass

    def disallow(self, stateA, stateB):
        pass

    def default_mover(self):
        pass


class MISTISNetwork(TISNetwork):
    def __init__(self):
        pass

    def default_mover(self):
        pass

