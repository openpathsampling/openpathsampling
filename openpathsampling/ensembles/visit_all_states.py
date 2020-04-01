import openpathsampling as paths
from openpathsampling.ensemble import EnsembleCache

def default_state_progress_report(n_steps, found_states, all_states,
                                  timestep=None):
    """
    Default progress reporter for VisitAllStatesEnsemble.

    Note that it is assumed that all states have been named.

    Parameters
    ----------
    n_steps : int
        number of MD frames generated so far
    found_states : iterable
        the set of states that have been found
    all_states : iterable
        the set of all states of interest
    timestep : float or quantity
        the timestep (optional). If given, the amount of time simulated will
        be reported along with the number of MD frames.

    Returns
    -------
    str :
        formatted string with information about progress so far
    """
    report_str = "Ran {n_steps} frames"
    if timestep is not None:
        report_str += " [{}]".format(str(n_steps * timestep))
    report_str += (". Found states [{found_states}]. "
                   "Looking for [{missing_states}].\n")
    found_states_str = ",".join([s.name for s in found_states])
    # list comprehension instead of sets (to preseve order)
    missing_states = [s for s in all_states if s not in found_states]
    missing_states_str = ",".join([s.name for s in missing_states])
    return report_str.format(n_steps=n_steps,
                             found_states=found_states_str,
                             missing_states=missing_states_str)


class VisitAllStatesEnsemble(paths.WrappedEnsemble):
    """
    Ensemble to create trajectories connected all states, giving progress.

    Parameters
    ----------
    states : list of :class:`.Volume`
        States this should visit
    progress : string or 2-tuple of callable
        How to report progress. This is used to define a method to create
        the report string and a method to emit the progress report.  Allowed
        string values are 'default', which uses
        :method:`.default_state_progress_report` to create the string and
        :method:`.refresh_output` to emit it, and 'silent', which will
        create the report as with the default, but not emit it.
        Using ``None`` will prevent the report from being created. You can
        also return a 2-tuple of callables, where the first formats the
        string (and must match the signature of
        :method:`.default_state_progress_report`) and the second tells how
        out emit the report.
    timestep : float or quantity
        Simulation timestep. See :method:`.default_state_progress_report`.

    Attributes
    ----------
    report_frequency : int
        how many frames between reports (default 10)
    progress_formatter : callable
        the method that creates the report string
    progress_emitter : callable
        the method that emits the report
    """
    def __init__(self, states, progress='default', timestep=None):
        self.states = states
        self.all_states = paths.join_volumes(states)
        all_states_ens = paths.join_ensembles([paths.AllOutXEnsemble(s)
                                               for s in states])
        ensemble = paths.SequentialEnsemble([
            all_states_ens,
            paths.AllInXEnsemble(self.all_states) & paths.LengthEnsemble(1)
        ])
        super(VisitAllStatesEnsemble, self).__init__(ensemble)
        self.timestep = timestep
        self.report_frequency = 10
        self.progress_formatter, self.progress_emitter = \
                self._progress_indicator(progress)
        self.cache = EnsembleCache(direction=+1)
        self._reset_cache_contents()

    def _reset_cache_contents(self):
        self.cache.contents['found_states'] = set([])

    @property
    def found_states(self):
        return self.cache.contents['found_states']

    @staticmethod
    def _progress_indicator(progress):
        """parse the input ``progress`` to get the actual objects"""
        indicator_dict = {
            None: (None, lambda x: None),
            'default': (default_state_progress_report,
                        paths.tools.refresh_output),
            'silent': (default_state_progress_report,
                       lambda x: None),
        }
        try:
            indicator = indicator_dict[progress]
        except KeyError:
            indicator = progress
        return indicator

    def _state_for_frame(self, snapshot):
        """get the state associated with a given frame"""
        in_states = [state for state in self.states if state(snapshot)]
        if len(in_states) > 1:
            raise RuntimeError("Frame is in more than one state.")
        elif len(in_states) == 1:
            state = set(in_states)
        else:
            state = set([])
        return state

    def progress_report(self, trajectory):
        """Convenience to get the progress report string for a trajectory.

        Parameters
        ----------
        trajectory : :class:`.Trajectory`
            input trajectory

        Returns
        -------
        str :
            progress report string
        """
        return self.progress_formatter(n_steps=len(trajectory) - 1,
                                       timestep=self.timestep,
                                       found_states=self.found_states,
                                       all_states=self.states)

    def _update_for_progress(self, trajectory, frame_number):
        """update cache of visited states; report if suitable"""
        len_traj = len(trajectory)
        frame = trajectory[frame_number]
        self.found_states.update(self._state_for_frame(trajectory[-1]))
        # minus 1 to account for the fact that the first step doesn't count
        if (len_traj - 1) % self.report_frequency == 0:
            report_string = self.progress_report(trajectory)
            self.progress_emitter(report_string)

    def can_append(self, trajectory, trusted=False):
        return_value = super(VisitAllStatesEnsemble, self).can_append(
            trajectory=trajectory,
            trusted=trusted
        )

        # TODO: can this be simplified and moved entirely into the cache?
        reset = self.cache.check(trajectory)
        if reset:
            self._reset_cache_contents()

        if self.progress_formatter:
            frames = [-1] if trusted else list(range(len(trajectory)))
            for frame in frames:
                self._update_for_progress(trajectory, frame_number=frame)

        if not return_value and self.progress_formatter:
            report_string = self.progress_report(trajectory)
            self.progress_emitter(report_string)

        return return_value

    strict_can_append = can_append

    def can_prepend(self, trajectory, trusted=False):
        raise NotImplementedError("prepend methods not implemented for ",
                                  "VisitAllStateEnsemble")

    strict_can_prepend = can_prepend

    def __call__(self, trajectory, candidate=False):
        if not candidate:
            self._reset_cache_contents()
            for frame in trajectory:
                self.found_states.update(self._state_for_frame(frame))

        return self.found_states == set(self.states)
