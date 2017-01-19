import openpathsampling as paths

from openpathsampling.netcdfplus import StorableNamedObject

class ChannelAnalysis(StorableNamedObject):
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
        super(ChannelAnalysis, self).__init__()
        self.channels = channels
        if steps is None:
            steps = []
        self.replica = replica

        self._treat_multiples = 'all'
        self._results = {c: [] for c in self.channels.keys() + [None]}
        self._analyze(steps)

    # separate this because I think much of the code might be generalized
    # later where step_num could be something else
    @staticmethod
    def _step_num(step):
        """Return ordinal number for the given input object.

        Abstracted so that other things might replace it.

        Parameters
        ----------
        step : :class:`.MCStep`
            the step

        Returns
        -------
        int :
            MC cycle number
        """
        return step.mccycle

    def _analyze(self, steps):
        """Primary analysis routine.

        Converts the input steps to an internal ._results dictionary of
        channel name to list of (start, end) tuples for when that channel is
        occupied.

        Parameters
        ----------
        steps : iterable of :class:`.MCStep`
            the steps to analyze
        """
        # for now, this assumes only one ensemble per channel
        # (would like that to change in the future)
        prev_traj = None
        last_start = {c: None for c in self._results}
        for step in steps:
            step_num = self._step_num(step)
            traj = step.active[self.replica].trajectory
            if prev_traj is None:
                prev_result = {c: len(self.channels[c].split(traj)) > 0
                               for c in self.channels}
                prev_result[None] = not any(prev_result.values())
                for c in last_start:
                    if prev_result[c] is True:
                        last_start[c] = step_num
            # re-use previous if the trajectory hasn't changed
            if traj is prev_traj:
                result = prev_result
            else:
                result = {c: len(self.channels[c].split(traj)) > 0
                          for c in self.channels}
                result[None] = not any(result.values())
                changed = [c for c in result if result[c] != prev_result[c]]
                for c in changed:
                    if result[c] is True:
                        # switched from False to True: entered this label
                        last_start[c] = step_num
                    else:
                        # switched from True to False: exited this label
                        finish = step_num
                        self._results[c] += [(last_start[c], finish)]
                        last_start[c] = None
            prev_traj = traj
            prev_result = result
        # finish off any extras
        next_step = step_num + 1 # again, this can be changed
        for c in self._results:
            if last_start[c] is not None:
                if len(self._results[c]) > 0:
                    # don't do double it if it's already there
                    if self._results[c][-1][1] != step_num:
                        self._results[c] += [(last_start[c], next_step)]
                    # note: is the else: of the above even possible?
                    # namely, do we need the if statement? should test that
                else:
                    self._results[c] += [(last_start[c], next_step)]

    @property
    def treat_multiples(self):
        return self._treat_multiples

    @treat_multiples.setter
    def treat_multiples(self, value):
        value = value.lower()
        if value not in ['all', 'newest', 'oldest', 'multiple']:
            raise ValueError("Invalid value for treat_multiples: " +
                             str(value))
        self._treat_multiples = value

    @staticmethod
    def _expand_results(results):
        expanded = [(domain[0], domain[1], channel)
                    for channel in results for domain in results[channel]]
        return sorted(expanded, key=lambda tup: tup[0])

    @staticmethod
    def _labels_by_step_newest(expanded_results):
        relabeled = []
        previous = expanded_results[0]
        for current in expanded_results[1:]:
            relabeled += [(previous[0], current[0], previous[2])]
            previous = current
        relabeled += [expanded_results[-1]]
        return relabeled

    @staticmethod
    def _labels_by_step_oldest(results):
        pass

    @staticmethod
    def _label_by_step_multiple(results):
        pass


    def labels_by_step(self):
        pass

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
