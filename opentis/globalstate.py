from wrapper import storable

@storable
class GlobalState(dict):
    """
    Notes
    =====
    I would suggest to add the timestamp to the GlobalState object and store it separately and not include
    the timestamp
    """

    now = None       # globsl current state
    time = 0         # globl current time

    def __init__(self, size, ensembles, time = None):
        super(GlobalState, self).__init__()
        self.size = size
        self.ensembles = ensembles
        self.old = None
        self.samples = []
        self.clear()
        if time is None:
            time = GlobalState.time + 1

        if time > GlobalState.time:
            GlobalState.time = time

        self.time = time
        for ensemble in ensembles:
            self[ensemble] = None

        GlobalState.now = self

    def move(self, samples):
        """
        Returns a new GlobalState object that takes the current instance and applies the sample in the given order as updates.

        samples : list of Sample()
            The list of `Sample` objects used to update the current GlobalState

        Returns
        =======
        GlobalState()

        """
        globalstate = GlobalState(self.size, self.ensembles)
        globalstate.time = self.time + 1
        globalstate.old = self
        globalstate.samples = samples
        for sample in samples:
            globalstate.ensemble_set[sample.ensemble] = sample.trajectory

        GlobalState.now = globalstate

        return globalstate