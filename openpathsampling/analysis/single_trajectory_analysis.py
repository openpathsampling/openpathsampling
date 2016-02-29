import openpathsampling as paths
import numpy as np


class SingleTrajectoryAnalysis(object):
    """Analyze a trajectory or set of trajectories for transition properties.

    Attributes
    ----------
    dt : float
        time step between frames
    continuous_frames : dict
        dictionary mapping state to a list of number of frames continuously
        in that state from the analyzed trajectories.
    lifetime_frames : dict
        dictionary mapping state to a list of the number of frames to the
        trajectory lengths for calculating the lifetime. See Notes for more.
    transition_frames : dict
        dictionary mapping the transition tuple (initial_state, final_state)
        to a list of the number of frames involves in the transition and not
        in either state.
    flux_frames : dict
        TODO
    continuous_times : dict
        As with continuous frames, but durations multiplied by self.dt
    lifetimes : dict
        As with lifetime_frames, but durations multiplied by self.dt
    transitions_durations : dict
        As with transition_frames, but durations multiplied by self.dt

    Notes
    -----
    LIFETIME: the
    """
    def __init__(self, transition, dt=None):
        self.transition = transition
        self.dt = dt
        self.stateA = transition.stateA
        self.stateB = transition.stateB
        self.continuous_frames = {self.stateA: np.array([]),
                                  self.stateB: np.array([])}
        self.lifetime_frames = {self.stateA: np.array([]),
                                self.stateB: np.array([])}
        self.transition_duration_frames = {
            (self.stateA, self.stateB): np.array([]),
            (self.stateB, self.stateA): np.array([])
        }

        self.flux_frames = {self.stateA: {}, self.stateB: {}}

    @property
    def continuous_times(self):
        if self.dt is None: # pragma: no-cover
            # TODO: this might become a logger.warn
            raise RuntimeError("No time delta set")
        return {k : self.continuous_frames[k]*self.dt 
                for k in self.continuous_frames.keys()}

    @property
    def lifetimes(self):
        if self.dt is None: # pragma: no-cover
            # TODO: this might become a logger.warn; use dt=1 otherwise
            raise RuntimeError("No time delta set")
        return {k : self.lifetime_frames[k]*self.dt 
                for k in self.lifetime_frames.keys()}

    @property
    def transition_duration(self):
        if self.dt is None: # pragma: no-cover
            # TODO: this might become a logger.warn; use dt=1 otherwise
            raise RuntimeError("No time delta set")
        return {k : self.transition_duration_frames[k]*self.dt 
                for k in self.transition_duration_frames.keys()}


    def analyze_continuous_time(self, trajectory, state):
        # convert to python list for append functionality
        continuous_frames = self.continuous_frames[state].tolist()
        ensemble = paths.AllInXEnsemble(state)
        segments = ensemble.split(trajectory)
        lengths = [len(seg) for seg in segments]
        continuous_frames += lengths
        # convert back to numpy to use as distribution
        self.continuous_frames[state] = np.array(continuous_frames)

    def analyze_lifetime(self, trajectory, state):
        # convert to python list for append functionality
        lifetime_frames = self.lifetime_frames[state].tolist()
        other_state = list(set([self.stateA, self.stateB]) - set([state]))[0]
	ensemble_BAB = paths.SequentialEnsemble([
	    paths.AllInXEnsemble(other_state) & paths.LengthEnsemble(1),
	    paths.PartInXEnsemble(state) & paths.AllOutXEnsemble(other_state),
	    paths.AllInXEnsemble(other_state) & paths.LengthEnsemble(1)
	])
	ensemble_AB = paths.SequentialEnsemble([
	    paths.AllInXEnsemble(state) & paths.LengthEnsemble(1),
	    paths.OptionalEnsemble(paths.AllOutXEnsemble(other_state)),
	    paths.AllInXEnsemble(other_state) & paths.LengthEnsemble(1)
	])
        BAB_split = ensemble_BAB.split(trajectory)
        AB_split = [ensemble_AB.split(part)[0] for part in BAB_split]
        lifetime_frames += [len(subtraj)-1 for subtraj in AB_split]

        # convert back to numpy to use as distribution
        self.lifetime_frames[state] = np.array(lifetime_frames)

    def analyze_transition_duration(self, trajectory, stateA, stateB):
        # we define the transitions ensemble just in case the transition is,
        # e.g., fixed path length TPS. We want flexible path length ensemble
        transition_ensemble = paths.SequentialEnsemble([
            paths.AllInXEnsemble(stateA) & paths.LengthEnsemble(1),
            paths.OptionalEnsemble( # optional to allow instantaneous hops
                paths.AllOutXEnsemble(stateA) & paths.AllOutXEnsemble(stateB)
            ),
            paths.AllInXEnsemble(stateB) & paths.LengthEnsemble(1)
        ])
        transition_segments = transition_ensemble.split(trajectory)
        self.transition_duration_frames[(stateA, stateB)] = np.append(
            self.transition_duration_frames[(stateA, stateB)],
            [len(seg)-2 for seg in transition_segments]
        )

    def analyze_flux(self, trajectory, state):
        pass

    def add_frames(self, trajectory):
        for state in [self.stateA, self.stateB]:
            self.analyze_continuous_time(trajectory, state)
            self.analyze_lifetime(trajectory, state)
            self.analyze_flux(trajectory, state)

    def analyze(self, trajectory):
        self.add_frames(trajectory)
        return self

    def summary(self):
        # return a nice text/pandas output
        pass
