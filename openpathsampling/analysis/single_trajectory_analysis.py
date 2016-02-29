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
        self.reset_analysis()

    def reset_analysis(self):
        self.continuous_segments = {self.stateA: [], self.stateB: []}
        self.lifetime_segments = {self.stateA: [], self.stateB: []}
        self.transition_segments = {(self.stateA, self.stateB): [], 
                                    (self.stateB, self.stateA): []}
        self.flux_segments = {self.stateA: {'in': [], 'out': []},
                              self.stateB: {'in': [], 'out': []}}

    @property
    def continuous_frames(self):
        return {k: np.array([len(seg) for seg in self.continuous_segments[k]])
                for k in self.continuous_segments.keys()}
    
    @property
    def continuous_times(self):
        if self.dt is None: # pragma: no-cover
            # TODO: this might become a logger.warn
            raise RuntimeError("No time delta set")
        continuous_frames = self.continuous_frames
        return {k : continuous_frames[k]*self.dt 
                for k in continuous_frames.keys()}

    @property
    def lifetime_frames(self):
        return {k: np.array([len(seg) for seg in self.lifetime_segments[k]])
                for k in self.lifetime_segments.keys()}
                             
    @property
    def lifetimes(self):
        if self.dt is None: # pragma: no-cover
            # TODO: this might become a logger.warn; use dt=1 otherwise
            raise RuntimeError("No time delta set")
        lifetime_frames = self.lifetime_frames
        return {k : lifetime_frames[k]*self.dt 
                for k in lifetime_frames.keys()}

    @property
    def transition_duration_frames(self):
        return {k: np.array([len(seg) for seg in self.transition_segments[k]])
                for k in self.transition_segments.keys()}

    @property
    def transition_duration(self):
        if self.dt is None: # pragma: no-cover
            # TODO: this might become a logger.warn; use dt=1 otherwise
            raise RuntimeError("No time delta set")
        transition_duration_frames = self.transition_duration_frames
        return {k : transition_duration_frames[k]*self.dt 
                for k in transition_duration_frames.keys()}


    def analyze_continuous_time(self, trajectory, state):
        ensemble = paths.AllInXEnsemble(state)
        self.continuous_segments[state] += ensemble.split(trajectory)

    def analyze_lifetime(self, trajectory, state):
        # convert to python list for append functionality
        #lifetime_frames = self.lifetime_frames[state].tolist()
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
        self.lifetime_segments[state] += [subtraj[0:-1] 
                                          for subtraj in AB_split]
        #lifetime_frames += [len(subtraj)-1 for subtraj in AB_split]

        # convert back to numpy to use as distribution
        #self.lifetime_frames[state] = np.array(lifetime_frames)

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
        self.transition_segments[(stateA, stateB)] += [
            seg[1:-1] for seg in transition_ensemble.split(trajectory)
        ]

    def analyze_flux(self, trajectory, state):
        pass

    def analyze(self, trajectories):
        # TODO: I hate using isinstance, but I don't see anoher way
        if isinstance(trajectories, paths.Trajectory):
            trajectories = [trajectories]
        for traj in trajectories:
            for state in [self.stateA, self.stateB]:
                self.analyze_continuous_time(traj, state)
                self.analyze_lifetime(traj, state)
                self.analyze_flux(traj, state)
        # return self so we can init and analyze in one line
        return self

    def summary(self):
        # return a nice text/pandas output
        pass
