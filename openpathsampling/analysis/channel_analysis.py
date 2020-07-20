import collections
import pandas as pd
import numpy as np

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
        self._results = {c: [] for c in list(self.channels.keys()) + [None]}
        if len(steps) > 0:
            self._analyze(steps)

    def to_dict(self):
        return {
            'results': self._results,
            'treat_multiples': self._treat_multiples,
            'channels': self.channels,
            'replica': self.replica
        }

    @classmethod
    def from_dict(cls, dct):
        obj = cls(steps=None,
                  channels=dct['channels'],
                  replica=dct['replica'])
        obj._results = dct['results']
        obj._treat_multiples = dct['treat_multiples']
        return obj

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
        """
        string : method for handling paths that match multiple channels

        Allowed values are:

        * 'newest': use the most recent channel entered
        * 'oldest': use the least recent channel entered
        * 'multiple': treat multiple channels as a new type of channel, e.g.,
          'a' and 'b' becomes 'a,b'
        * 'all': treat each channel individually, despite overlaps. For
          switching, this is the same as ???. For status, this is the
          same as 'multiple'

        """
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
        """
        Takes ._results dict and makes it into chronological list of events

        Note
        ----
            The output of this is in terms of "channel events." It doesn't
            do anything to ensure that only one channel is defined at a
            given time --- so subsequent events can include overlapping
            times. Other functions relabel by time.

        See also
        --------
            _labels_by_step_newest
            _labels_by_step_oldest
            _labels_by_step_multiple

        Parameters
        ----------
        results : dict of {str: [(int, int), ...]}
            the results dictionary. The keys are the channel names, and the
            values are a list of tuples representing the start and finish
            step numbers for the range of steps while in this channel

        Returns
        -------
        list of 3-tuples (int, int, frozenset) :
            the "events": each event is the tuple of start step, finish
            step, and channel name (as a frozenset containing one string),
            sorted according to the start step.
        """
        expanded = [(domain[0], domain[1], frozenset([channel]))
                    for channel in results for domain in results[channel]]
        return sorted(expanded, key=lambda tup: tup[0])

    @staticmethod
    def _labels_by_step_newest(expanded_results):
        """
        Makes one channel per step, based on most recent channel entered.

        See also
        --------
            _expand_results
            _labels_by_step_oldest
            _labels_by_step_multiple
            labels_by_step

        Parameters
        ----------
        expanded_results : list of 3-tuples (int, int, frozenset)
            input events; the output of _expand_results (see details there)

        Returns
        -------
        list of 3-tuples (int, int, frozenset)
            events with ranges such that there is only one channel at any
            given step number, and that channel is the most recent entered
        """
        relabeled = []
        previous = expanded_results[0]
        for current in expanded_results[1:]:
            relabeled += [(previous[0], current[0], previous[2])]
            previous = current
        relabeled += [expanded_results[-1]]
        return relabeled

    @staticmethod
    def _labels_by_step_oldest(expanded_results):
        """
        Makes one channel per step, based on least recent channel entered.

        See also
        --------
            _expand_results
            _labels_by_step_newest
            _labels_by_step_multiple
            labels_by_step

        Parameters
        ----------
        expanded_results : list of 3-tuples (int, int, frozenset)
            input events; the output of _expand_results (see details there)

        Returns
        -------
        list of 3-tuples (int, int, frozenset) :
            events with ranges such that there is only one channel at any
            given step number, and that channel is the least recent entered
        """
        relabeled = []
        previous = expanded_results[0]
        for current in expanded_results[1:]:
            if current[1] > previous[1]:
                # ends after last one ended
                # if this isn't true, this one gets skipped
                # if it is true, then previous is used
                relabeled += [previous]
                # save the new starting point
                previous = (previous[1], current[1], current[2])
            # NOTE: Tests include a case for the implicit "else" here, but
            # for some reason an empty `else: pass` doesn't show up as
            # covering the `pass` line; so removed them
        if relabeled[-1] != previous:
            relabeled += [previous]
        else:
            pass # for testing
        return relabeled

    @staticmethod
    def _labels_by_step_multiple(expanded_results):
        """Makes one channel label per step, combining all active channels.

        See also
        --------
            _expand_results
            _labels_by_step_newest
            _labels_by_step_oldest
            labels_by_step

        Parameters
        ----------
        expanded_results : list of 3-tuples (int, int, frozenset)
            input events; the output of _expand_results (see details there)

        Returns
        -------
        list of 3-tuples (int, int, frozenset) :
            events such that there is only one event at any given step
            number, with the channel label as the set of all active channels
        """
        relabeled = []
        # start events are times when a channel is added to the active
        # finish events are when channel is removed from the active
        # both are dicts of time to a set of channels
        start_events = collections.defaultdict(set)
        finish_events = collections.defaultdict(set)
        for event in expanded_results:
            start_events[event[0]] |= set(event[2])
            finish_events[event[1]] |= set(event[2])

        all_event_steps = set(start_events.keys()) | set(finish_events.keys())
        active_channels = set([])
        prev_step_num = None
        # note to self: this is some elegant freaking code
        for step_num in sorted(list(all_event_steps)):
            if prev_step_num is not None:
                relabeled += [(prev_step_num, step_num,
                               frozenset(active_channels))]

            # defaultdict gives empty if doesn't exist
            active_channels -= finish_events[step_num]
            active_channels |= start_events[step_num]

            prev_step_num = step_num

        return relabeled

    def labels_by_step(self, treat_multiples=None):
        """
        Prepare internally stored results for primary analysis routines.

        Note
        ----
            The results of this depend on the value of ``treat_multiples``.
            In fact, this method is just a switch for the specific
            ``treat_multiples`` values.

        See also
        --------
            _labels_by_step_newest
            _labels_by_step_oldest
            _labels_by_step_multiple

        Parameters
        ----------
            treat_multiples : 'all', 'newest', 'oldest',  'multiple', or None
                method to prepare output; see documentation on
                ``treat_multiples`` for details. Default is `None`, which
                uses the value in ``self.treat_multiples``.

        Returns
        -------
        list of 3-tuples (int, int, frozenset) :
            events such that there is only one event at any give step
        """
        if treat_multiples is None:
            treat_multiples = self.treat_multiples
        expanded_results = self._expand_results(self._results)
        method = {
            'all': lambda x: x,
            'newest': self._labels_by_step_newest,
            'oldest': self._labels_by_step_oldest,
            'multiple': self._labels_by_step_multiple
        }[treat_multiples]
        return method(expanded_results)

    @staticmethod
    def _labels_as_sets_sort_function(label):
        """Sort function for labels.

        The input labels are frozensets of lists of strings. The sort order
        is first by number of items in the list, and then by the list items.
        A list of None will be sorted into the first place.

        Parameters
        ----------
        label: frozenset of list of (string or None)
            input label

        Returns
        -------
        list:
            first element is the length of the input set, followed by
            the input as a sorted list
        """
        label_list = list(label)
        if None in label_list:
            has_None = [None]
            label_list.remove(None)
        else:
            has_None = []
        ll = sorted(label_list)
        return [len(ll)] + has_None + ll

    @staticmethod
    def label_to_string(label):
        """Convert set of string/None to comma-separated string.

        For example, frozenset(['c', 'a']) becomes 'a,c' (no space).

        Parameters
        ----------
        label: frozenset of list of (string or None)
            input label

        Returns
        -------
        string:
            the string for this label
        """
        # separated for reusability
        return ",".join(sorted([str(l) for l in list(label)]))

    @property
    def switching_matrix(self):
        """
        pandas.DataFrame :
            number of switches from one channel to another. Depends on
            ``treat_multiples``, see details there.
        """
        labeled_results = self.labels_by_step()
        labels_in_order = [ll[2] for ll in labeled_results]
        labels_set = set(labels_in_order)
        sorted_set_labels = sorted(list(labels_set),
                                   key=self._labels_as_sets_sort_function)
        sorted_labels = [self.label_to_string(e) for e in sorted_set_labels]
        switches = [(self.label_to_string(labels_in_order[i]),
                     self.label_to_string(labels_in_order[i+1]))
                    for i in range(len(labeled_results)-1)]
        switch_count = collections.Counter(switches)
        df = pd.DataFrame(index=sorted_labels, columns=sorted_labels)
        for switch in switch_count:
            df.at[switch[0], switch[1]] = switch_count[switch]

        df = df.fillna(0)
        return df

    @property
    def residence_times(self):
        """
        Dict[string, List[int]] :
            number of steps spent in each channel for each "stay" in that
            channel; allows calculations of distribution properties. Depends
            on ``treat_multiples``, see details there.
        """
        labeled_results = self.labels_by_step()
        durations = [(self.label_to_string(step[2]), step[1] - step[0])
                     for step in labeled_results]
        results = collections.defaultdict(list)
        for dur in durations:
            results[dur[0]] += [dur[1]]
        return results

    @property
    def total_time(self):
        """
        Dict[string, int] :
            total number of steps spent in each channel for each "stay" in
            that channel. Depends on ``treat_multiples``, see details there.
        """
        residences = self.residence_times
        results = collections.defaultdict(int)
        for channel in residences:
            results[channel] = sum(residences[channel])
        return results

    def status(self, step_number):
        """Reports which channel(s) are associated with a given step number.

        Note
        ----
            Results will depend on the value of ``treat_multiples``. See
            details in the documentation for that.

        Parameters
        ----------
        step_number : int
            the step number of interest

        Returns
        -------
        string :
            the string label for the channel(s)
        """
        treat_multiples = self.treat_multiples
        if self.treat_multiples == 'all':
            treat_multiples = 'multiple'
        labeled_results = self.labels_by_step(treat_multiples)
        for step in labeled_results:
            if step[0] <= step_number < step[1]:
                return self.label_to_string(step[2])
        raise RuntimeError("Step " + str(step_number) + " outside of range."
                           + " Max step: " + str(labeled_results[-1][1]))
