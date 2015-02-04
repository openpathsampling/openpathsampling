class Transition(object):
    def __init__(self, stateA, stateB, interfaces, storage=None):
        pass

    @property
    def sample_dict(self):
        pass

    def load_samples(self, storage=None):
        pass


class CrossingProbabilities(object):
    def __init__(self, transition, order_parameter, combine="wham"):
        pass


class MinusMoveFlux(object):
    def __init__(self, storage, state, through_interfaces=[]):
        pass

    @property
    def sample_dict(self):
        pass

    def load_samples(self, storage=None):
        pass


class TransitionRate(object):
    def __init__(self, transition):
        pass


class ReplicaFlow(object):
    def __init__(self, top, bottom):
        pass

    def __call__(self, storage):
        pass


class PathEnsembleStats(object):
    def __init__(self, ensemble):
        pass


class MultipleSetMinusSwitching(object):
    def __init__(self, state, storage=None):
        pass



