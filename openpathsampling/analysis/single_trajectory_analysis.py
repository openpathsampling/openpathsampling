import openpathsampling as paths
import numpy as np

class TrajectorySegmentContainer(object):
    def __init__(self, segments, dt=None):
        self._segments = segments
        self.dt = dt

    @property
    def n_frames(self):
        return np.array([len(seg) for seg in self._segments])

    @property
    def times(self):
        if self.dt is None:
            raise RuntimeError("No time delta set")
            # TODO: this might become a logger.warn
        return np.array([len(seg)*self.dt for seg in self._segments])

    def __add__(self, other):
        if self.dt != other.dt:
            raise RuntimeError(
                "Different time steps in TrajectorySegmentContainers."
            )
        return TrajectorySegmentContainer(self._segments + other._segments,
                                          self.dt)

    def __iadd__(self, other):
        # in this case, we ignore dt
        self._segments += other._segments
        return self

    def __len__(self):
        return len(self._segments)

    def __contains__(self, item):
        return item in self._segments

    def __iter__(self):
        return self._segments.__iter__()

    def __getitem__(self, key):
        return self._segments[key]

    def __setitem__(self, key, value):
        raise TypeError("TrajectorySegmentContainer is immutable")

    # intentionally do not support __mul__ & related
    

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
        dictionary mapping each state to a dictionary of with keys 'in' (for
        frames in the state) and 'out' (for frames outside the state)
    continuous_times : dict
        As with continuous frames, but durations multiplied by self.dt
    lifetimes : dict
        As with lifetime_frames, but durations multiplied by self.dt
    transitions_durations : dict
        As with transition_frames, but durations multiplied by self.dt

    """
    def __init__(self, transition, dt=None):
        self.transition = transition
        self.dt = dt
        self.stateA = transition.stateA
        self.stateB = transition.stateB
        self.reset_analysis()

    def reset_analysis(self):
        """Reset the analysis by emptying all saved segments."""
        stateA = self.stateA
        stateB = self.stateB
        dt = self.dt
        self.continuous_segments = {
            stateA: TrajectorySegmentContainer(segments=[], dt=dt), 
            stateB: TrajectorySegmentContainer(segments=[], dt=dt)
        }
        self.lifetime_segments = {
            stateA: TrajectorySegmentContainer(segments=[], dt=dt), 
            stateB: TrajectorySegmentContainer(segments=[], dt=dt)
        }
        self.transition_segments = {
            (stateA, stateB): TrajectorySegmentContainer(segments=[], dt=dt),
            (stateB, stateA): TrajectorySegmentContainer(segments=[], dt=dt)
        }
        self.flux_segments = {
            stateA: {'in': TrajectorySegmentContainer(segments=[], dt=dt), 
                     'out': TrajectorySegmentContainer(segments=[], dt=dt)},
            stateB: {'in': TrajectorySegmentContainer(segments=[], dt=dt), 
                     'out': TrajectorySegmentContainer(segments=[], dt=dt)}
        }

    @property
    def continuous_frames(self):
        return {k: np.array([len(seg) for seg in self.continuous_segments[k]])
                for k in self.continuous_segments.keys()}
    
    @property
    def continuous_times(self):
        if self.dt is None: # pragma: no cover
            raise RuntimeError("No time delta set")
            # TODO: this might become a logger.warn
        continuous_frames = self.continuous_frames
        return {k : continuous_frames[k]*self.dt 
                for k in continuous_frames.keys()}

    @property
    def lifetime_frames(self):
        return {k: np.array([len(seg) for seg in self.lifetime_segments[k]])
                for k in self.lifetime_segments.keys()}
                             
    @property
    def lifetimes(self):
        if self.dt is None: # pragma: no cover
            raise RuntimeError("No time delta set")
            # TODO: this might become a logger.warn; use dt=1 otherwise
        lifetime_frames = self.lifetime_frames
        return {k : lifetime_frames[k]*self.dt 
                for k in lifetime_frames.keys()}

    @property
    def transition_duration_frames(self):
        return {k: np.array([len(seg) for seg in self.transition_segments[k]])
                for k in self.transition_segments.keys()}

    @property
    def transition_duration(self):
        if self.dt is None: # pragma: no cover
            raise RuntimeError("No time delta set")
            # TODO: this might become a logger.warn; use dt=1 otherwise
        transition_duration_frames = self.transition_duration_frames
        return {k : transition_duration_frames[k]*self.dt 
                for k in transition_duration_frames.keys()}


    def analyze_continuous_time(self, trajectory, state):
        """Analysis to obtain continuous times for given state.

        Parameters
        ----------
        trajectory : :class:`.Trajectory`
            trajectory to analyze
        state : :class:`.Volume`
            state volume to characterize. Must be one of the states in the
            transition
        """
        ensemble = paths.AllInXEnsemble(state)
        segments = ensemble.split(trajectory, overlap=0)
        new_container = TrajectorySegmentContainer(segments, self.dt)
        self.continuous_segments[state] += new_container
    
    @staticmethod
    def get_lifetime_segments(trajectory, from_vol, to_vol, forbidden=None,
                              padding=[0, -1]):
        """General script to get lifetimes.

        Lifetimes for a transition between volumes are used in several other
        calculations: obviously, the state lifetime, but also the flux
        through an interface. This is a generic function to calculate that.

        Parameters
        ----------
        trajectory : :class:`.Trajectory`
            trajectory to analyze
        from_vol : :class:`.Volume`
            the volume for which this represents the lifetime: the
            trajectory segments returned are associated with the lifetime of
            `from_vol`
        to_vol : :class:`.Volume`
            the volume which indicates the end of the lifetime: a frame in
            this volume means the trajectory is no longer associated with
            `from_vol`
        forbidden : :class:`.Volume`
            if a frame is in `forbidden`, it cannot be part of the lifetime
            of `from_vol`. This isn't needed in 2-state lifetime
            calculations; however, it is useful to exclude other states
            from a flux calculation
        padding : list
            adjusts which frames are returned as list indices. That is, the
            returned segments are `full_segment[padding[0]:padding[1]]`.
            The `full_segment`s are the segments from (and including) each
            first frame in `from_vol` (after a visit to `to_vol`) until (and
            including) the first frame in `to_vol`. To get the full segment
            as output, use `padding=[None, None]`. The default is to remove
            the final frame (`padding=[0, -1]`) so that it doesn't include
            the frame in `to_vol`.

        Returns
        -------
        :class:`.TrajectorySegmentContainer`
            the frames from (and including) each first entry from `to_vol`
            into `from_vol` until (and including) the next entry into
            `to_vol`, with no frames in `forbidden`, and with frames removed
            from the ends according to `padding`
        """
        if forbidden is None:
            forbidden = paths.EmptyVolume()
        ensemble_BAB = paths.SequentialEnsemble([
            paths.AllInXEnsemble(to_vol) & paths.LengthEnsemble(1),
            paths.PartInXEnsemble(from_vol) & paths.AllOutXEnsemble(to_vol),
            paths.AllInXEnsemble(to_vol) & paths.LengthEnsemble(1)
        ]) & paths.AllOutXEnsemble(forbidden)
        ensemble_AB = paths.SequentialEnsemble([
            paths.AllInXEnsemble(from_vol) & paths.LengthEnsemble(1),
            paths.OptionalEnsemble(paths.AllOutXEnsemble(to_vol)),
            paths.AllInXEnsemble(to_vol) & paths.LengthEnsemble(1)
        ])
        BAB_split = ensemble_BAB.split(trajectory)
        AB_split = [ensemble_AB.split(part)[0] for part in BAB_split]
        return [subtraj[padding[0]:padding[1]] for subtraj in AB_split]


    def analyze_lifetime(self, trajectory, state):
        """Analysis to obtain  lifetimes for given state.

        Parameters
        ----------
        trajectory : :class:`.Trajectory`
            trajectory to analyze
        state : :class:`.Volume`
            state volume to characterize. Must be one of the states in the
            transition
        """
        other_state = list(set([self.stateA, self.stateB]) - set([state]))[0]
        segments = self.get_lifetime_segments(
            trajectory=trajectory,
            from_vol=state,
            to_vol=other_state
        )
        self.lifetime_segments[state] += TrajectorySegmentContainer(segments,
                                                                    self.dt)

    def analyze_transition_duration(self, trajectory, stateA, stateB):
        """Analysis to obtain transition durations for given state.

        Parameters
        ----------
        trajectory : :class:`.Trajectory`
            trajectory to analyze
        stateA : :class:`.Volume`
            initial state volume for the transition
        stateB : :class:`.Volume`
            final state volume for the transition
        """
        # we define the transitions ensemble just in case the transition is,
        # e.g., fixed path length TPS. We want flexible path length ensemble
        transition_ensemble = paths.SequentialEnsemble([
            paths.AllInXEnsemble(stateA) & paths.LengthEnsemble(1),
            paths.OptionalEnsemble( # optional to allow instantaneous hops
                paths.AllOutXEnsemble(stateA) & paths.AllOutXEnsemble(stateB)
            ),
            paths.AllInXEnsemble(stateB) & paths.LengthEnsemble(1)
        ])
        segments = [seg[1:-1] for seg in transition_ensemble.split(trajectory)]
        new_container = TrajectorySegmentContainer(segments, self.dt)
        self.transition_segments[(stateA, stateB)] += new_container

    def analyze_flux(self, trajectory, state, interface=None):
        """Analysis to obtain flux segments for given state.

        Parameters
        ----------
        trajectory : :class:`.Trajectory`
            trajectory to analyze
        state : :class:`.Volume`
            state volume to characterize. Must be one of the states in the
            transition
        interface : :class:`.Volume` or None
            interface to calculate the flux through. If `None`, same as
            `state`
        """
        other = list(set([self.stateA, self.stateB]) - set([state]))[0]
        if interface is None:
            interface = state
        out_segments = self.get_lifetime_segments(
            trajectory=trajectory,
            from_vol=~interface,
            to_vol=state,
            forbidden=other,
            padding=[None, -1]
        )
        out_container = TrajectorySegmentContainer(out_segments, self.dt)
        self.flux_segments[state]['out'] += out_container
        in_segments = self.get_lifetime_segments(
            trajectory=trajectory,
            from_vol=state,
            to_vol=~interface,
            forbidden=other,
            padding=[None, -1]
        )
        in_container = TrajectorySegmentContainer(in_segments, self.dt)
        self.flux_segments[state]['in'] += in_container


    def analyze(self, trajectories):
        """Full analysis of a trajectory or trajectories.

        Parameters
        ----------
        trajectories : :class:`.Trajectory` or list of :class:`.Trajectory`
        """
        # TODO: I hate using isinstance, but I don't see another way
        if isinstance(trajectories, paths.Trajectory):
            trajectories = [trajectories]
        for traj in trajectories:
            for state in [self.stateA, self.stateB]:
                self.analyze_continuous_time(traj, state)
                self.analyze_lifetime(traj, state)
                self.analyze_flux(traj, state)
            self.analyze_transition_duration(traj, self.stateA, self.stateB)
            self.analyze_transition_duration(traj, self.stateB, self.stateA)
        # return self so we can init and analyze in one line
        return self

    # TODO: add a `summary` function to output a nice pandas frame or
    # something
