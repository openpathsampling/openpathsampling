import logging
import numpy as np
import pandas as pd

import openpathsampling as paths
from .path_simulator import PathSimulator, MCStep

logger = logging.getLogger(__name__)
init_log = logging.getLogger('openpathsampling.initialization')

# python 3 support
try:
    xrange
except NameError:
    xrange = range

class DirectSimulation(PathSimulator):
    """
    Direct simulation to calculate rates and fluxes.

    In practice, this is primarily used to calculate the flux if you want to
    do so without saving the entire trajectory. However, it will also save
    the trajectory, if you want it to.

    Parameters
    ----------
    storage : :class:`.Storage`
        file to store the trajectory in. Default is None, meaning that the
        trajectory isn't stored (also faster)
    engine : :class:`.DynamicsEngine`
        the engine for the molecular dynamics
    states : list of :class:`.Volume`
        states to look for transitions between
    flux_pairs : list of 2-tuples of ``(state, interface)``
        fluxes will calculate the flux out of `state` and through
        `interface` for each pair in this list
    initial_snapshot : :class:`.Snapshot`
        initial snapshot for the MD

    Attributes
    ----------
    transitions : dict with keys 2-tuple of paths.Volume, values list of int
        for each pair of states (from_state, to_state) as a key, gives the
        number of frames for each transition from the entry into from_state
        to entry into to_state
    rate_matrix : pd.DataFrame
        calculates the rate matrix, in units of per-frames
    fluxes : dict with keys 2-tuple of paths.Volume, values float
        flux out of state and through interface for each (state, interface)
        key pair
    n_transitions : dict with keys 2-tuple of paths.Volume, values int
        number of transition events for each pair of states
    n_flux_events : dict with keys 2-tuple of paths.Volume, values int
        number of flux events for each (state, interface) pair
    """
    def __init__(self, storage=None, engine=None, states=None,
                 flux_pairs=None, initial_snapshot=None):
        super(DirectSimulation, self).__init__(storage)
        self.engine = engine
        self.states = states
        self.flux_pairs = flux_pairs
        if flux_pairs is None:
            self.flux_pairs = []
        self.initial_snapshot = initial_snapshot
        self.save_every = 1

        # TODO: might set these elsewhere for reloading purposes?
        self.transition_count = []
        self.flux_events = {pair: [] for pair in self.flux_pairs}

    @property
    def results(self):
        return {'transition_count': self.transition_count,
                'flux_events': self.flux_events}

    def load_results(self, results):
        self.transition_count = results['transition_count']
        self.flux_events = results['flux_events']

    def run(self, n_steps):
        most_recent_state = None
        first_interface_exit = {p: -1 for p in self.flux_pairs}
        last_state_visit = {s: -1 for s in self.states}
        was_in_interface = {p: None for p in self.flux_pairs}
        local_traj = paths.Trajectory([self.initial_snapshot])
        self.engine.current_snapshot = self.initial_snapshot
        self.engine.start()
        for step in xrange(n_steps):
            frame = self.engine.generate_next_frame()

            # update the most recent state if we're in a state
            state = None  # no state at all
            for s in self.states:
                if s(frame):
                    state = s
            if state:
                last_state_visit[state] = step
                if state is not most_recent_state:
                    # we've made a transition: on the first entrance into
                    # this state, we reset the last_interface_exit
                    state_flux_pairs = [p for p in self.flux_pairs
                                        if p[0] == state]
                    for p in state_flux_pairs:
                        first_interface_exit[p] = -1
                    # if this isn't the first change of state, we add the
                    # transition
                    if most_recent_state:
                        self.transition_count.append((state, step))
                    most_recent_state = state

            # update whether we've left any interface
            for p in self.flux_pairs:
                state = p[0]
                interface = p[1]
                is_in_interface = interface(frame)
                # by line: (1) this is a crossing; (2) the most recent state
                # is correct; (3) this is the FIRST crossing
                first_exit_condition = (
                    not is_in_interface and was_in_interface[p]  # crossing
                    and state is most_recent_state  # correct recent state
                    and first_interface_exit[p] < last_state_visit[state]
                )
                if first_exit_condition:
                    first_exit = first_interface_exit[p]
                    # successful exit
                    if 0 < first_exit < last_state_visit[state]:
                        flux_time_range = (step, first_exit)
                        self.flux_events[p].append(flux_time_range)
                    first_interface_exit[p] = step
                was_in_interface[p] = is_in_interface

            if self.storage is not None:
                local_traj += [frame]

        self.engine.stop(local_traj)

        if self.storage is not None:
            self.storage.save(local_traj)

    @property
    def transitions(self):
        prev_state = None
        prev_time = None
        results = {}
        for (new_state, time) in self.transition_count:
            if prev_state is not None and prev_time is not None:
                lag = time - prev_time
                try:
                    results[(prev_state, new_state)] += [lag]
                except KeyError:
                    results[(prev_state, new_state)] = [lag]
            prev_state = new_state
            prev_time = time
        return results

    @property
    def rate_matrix(self):
        transitions = self.transitions
        try:
            time_per_step = self.engine.snapshot_timestep
        except AttributeError:
            time_per_step = 1.0
        total_time = {s: sum(sum((transitions[t] for t in transitions
                                  if t[0] == s), [])) * time_per_step
                      for s in self.states}

        rates = {t : len(transitions[t]) / total_time[t[0]]
                 for t in transitions}
        # rates = {t : 1.0 / np.array(transitions[t]).mean()
                 # for t in transitions}

        state_names = [s.name for s in self.states]
        rate_matrix = pd.DataFrame(columns=state_names, index=state_names)
        for t in rates:
            rate_matrix.at[t[0].name, t[1].name] = rates[t]
        return rate_matrix

    @property
    def fluxes(self):
        results = {}
        try:
            time_per_step = self.engine.snapshot_timestep
        except AttributeError:
            time_per_step = 1.0

        for p in self.flux_events:
            lags = [t[0] - t[1] for t in self.flux_events[p]]
            results[p] = 1.0 / np.mean(lags) / time_per_step
        return results

        # return {p : 1.0 / np.array(self.flux_events[p]).mean()
                # for p in self.flux_events}

    @property
    def n_transitions(self):
        transitions = self.transitions
        return {t : len(transitions[t]) for t in transitions}

    @property
    def n_flux_events(self):
        return {p : len(self.flux_events[p]) for p in self.flux_events}
