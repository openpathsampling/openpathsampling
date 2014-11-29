from wrapper import storable

@storable
class GlobalState(dict):
    """
    Notes
    =====
    I would suggest to add the timestamp to the GlobalState object and store
    it separately and not include the timestamp in the samples since these
    might be saved before a timestamp is set and thus we have mutable
    objects the might screw up the storage
    """

    current = None       # global current state
    now = 0         # global current time

    def __init__(self, ensembles, time=None):
        super(GlobalState, self).__init__()
        self._ensembles = ensembles
        self.old = None
        self.samples = []
        self.clear()
        if time is None:
            time = GlobalState.now + 1

        if time > GlobalState.now:
            GlobalState.now = time

        self.step = time
        for ensemble in self._ensembles:
            dict.__setitem__(self, ensemble, None)

        GlobalState.current = self

    def __setitem__(self, key, value):
        if key in self:
            if self[key] is None:
                dict.__setitem__(self, key, value)
            else:
                print 'Cannot set a trajectory for this Ensemble. Already set. ' \
                      'Trajectories can only be set once. Otherwise moves are necessary!'

        else:
            print 'Cannot change size of GlobalState. Ensemble are fixed at creation!'

    def __getitem__(self, item):
        if type(item) is int:
            return dict.__getitem__(self, self.ensembles[item])
        else:
            return dict.__getitem__(self, item)

    def __repr__(self):
        return repr([ self[ens] for ens in self._ensembles ])


    @property
    def size(self):
        """
        Returns the number of unique ensembles

        Returns
        =======
        int
            The number of ensembles / trajectories in the GlobalState
        """
        return len(self._ensembles)

    @property
    def ensembles(self):
        """
        Return an orderes list of ensembles in the GlobalState

        Returns
        =======
        list of Ensemble()
            The initial list of ensembles in the given order in the GlobalState
        """
        return self._ensembles

    @property
    def trajectories(self):
        """
        Return an ordered list of trajectories in the GlobalState()

        Returns
        =======
        list of Trajectory()
            The list of trajectories
        """
        return [ self[ensemble] for ensemble in self._ensembles]

    def apply_samples(self, samples):
        """
        Returns a new GlobalState object that takes the current instance and
        applies the samples in the given order as updates.  The samples will
        get the actual timestamp for later analysis.

        Parameters
        ==========
        samples : list of Sample()
            The list of `Sample` objects used to update the current GlobalState

        Returns
        =======
        GlobalState()

        """
        globalstate = GlobalState(self.ensembles)
        globalstate.now = self.step + 1
        globalstate.old = self
        self.samples = samples
        for ensemble in self.ensembles:
            globalstate[ensemble] = self[ensemble]

        for sample in samples:
            sample.time = self.step

            dict.__setitem__(globalstate, sample.ensemble, sample.trajectory)

        GlobalState.current = globalstate
        return globalstate

    def save_samples(self, storage):
        """
        Save all samples in the current GlobalState object. This should be
        called after a move has generated a new object since then all
        samples will get a timestamp that is associated with this

        Parameters
        ==========
        storage : Storage()
            the underlying netcdf file to be used for storage
        """
        map(storage.sample.save, self.samples)


class GlobalStateMover(object):
    """
    A GlobalStateMover will take a GlobalState object and create in a
    specific way a new GlobalState by generating the necessary samples,
    applying them and saving everything correctly.
    """
    def __init__(self, movers)
        pass

    def move(self, globalstate):
        pass
    

class GlobalMixedMover(GlobalStateMover):
    def __init__(self, movers, probabilities):
        pass

    def move(self, globalstate):
        pass


class GeneralPathMover(GlobalStateMover):
    # TODO: takes a list of list of moves; default is that each individual
    # move has the same probability ... 
    # list of moves becomes GlobalMover, list of GlobalMovers becomes
    # GlobalMixedMover
    pass
