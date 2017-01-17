import openpathsampling as paths

class ChannelAnalysis(paths.StorableNamedObject):
    """Analyze path sampling simulation for multiple channels.

    User defines several channels (e.g., mechanisms) as :class:`.Ensemble`
    objects. This checks which channels each path satisfies, and provides
    analysis of switching and residence.

    Parameters
    ----------
    steps : iterable of :class:`.MCStep`
        the steps to analyze
    channels: dict of {string: :class:`.Ensemble`}
        names (keys) and ensembles (values) representing subtrajectories of the
        channels of interest
    replica: int
        replica ID to analyze from the steps, default is 0.
    """
    def __init__(self, steps, channels, replica=0):
        self.channels = channels
        if steps is None:
            steps = []
        self.replica = replica

        self.treat_multiples = 'unique'
        self._results = {c: [] for c in self.channels}
        self._analyze(steps)

    # separate this because I think much of the code might be generalized
    # later where step_num could be something else
    @staticmethod
    def _step_num(step):
        return step.mccycle

    def _analyze(steps):
        """Primary analysis routine.

        Parameters
        ----------
        steps : iterable of :class:`.MCStep`
            the steps to analyze
        """
        # for now, this assumes only one ensemble per channel
        # (would like that to change in the future)
        prev_traj = None
        prev_result = None
        last_start = {c: None for c in self.channels}
        for step in steps:
            traj = step.active[replica].trajectory
            # re-use previous if the trajectory hasn't changed
            if traj is prev_traj:
                result = prev_result
            else:
                result = {c: len(self.channels[c].split(traj)) > 0
                          for c in self.channels}
                changed = [c for c in result if result[c] != prev_result[c]]
                for c in changed:
                    if result[c] is True:
                        # switched from False to True: entered this label
                        last_start[c] = self._step_num(step)
                    else:
                        # switched from True to False: exited this label
                        finish = self._step_num(step)
                        self._results[c] += (last_start[c], finish)
                        last_start[c] = None

    @property
    def treat_multiples(self):
        return self._treat_multiples

    @treat_multiples.setter
    def treat_multiples(self, value):
        self._treat_multiples = value

    @property
    def switching_matrix(self):
        pass

    @property
    def residence_times(self):
        pass

    @property
    def total_time(self):
        pass

    def status(self, step_number):
        """Which channels were active at a given step number"""
        pass
