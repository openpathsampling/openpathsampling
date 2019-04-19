import openpathsampling as paths
from openpathsampling.ensemble import EnsembleCache

def default_state_progress_report(n_steps, found_states, all_states,
                                  timestep=None):
    report_str = "Ran {n_steps} steps"
    if timestep is not None:
        report_str += " [{}]".format(str(n_steps * timestep))
    report_str += (". Found states [{found_states}]. "
                   "Looking for [{missing_states}].")
    found_states_str = ",".join([s.name for s in found_states])
    # list comprehension instead of sets (to preseve order)
    missing_states = [s for s in all_states if s not in found_states]
    missing_states_str = ",".join([s.name for s in missing_states])
    return report_str.format(n_steps=n_steps,
                             found_states=found_states_str,
                             missing_states=missing_states_str)


class VisitAllStatesEnsemble(paths.WrappedEnsemble):
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
        in_states = [state for state in self.states if state(snapshot)]
        if len(in_states) > 1:
            raise RuntimeError("Frame is in more than one state.")
        elif len(in_states) == 1:
            state = set(in_states)
        else:
            state = set([])
        return state

    def progress_report(self, trajectory):
        return self.progress_formatter(n_steps=len(trajectory) - 1,
                                       timestep=self.timestep,
                                       found_states=self.found_states,
                                       all_states=self.states)

    def _update_for_progress(self, trajectory, frame_number):
        len_traj = len(trajectory)
        frame = trajectory[frame_number]
        self.found_states.update(self._state_for_frame(trajectory[-1]))
        if len_traj - 1 % self.report_frequency == 0:
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
