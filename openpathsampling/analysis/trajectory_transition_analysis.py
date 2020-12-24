import openpathsampling as paths
import numpy as np

class TrajectorySegmentContainer(object):
    """Container object to analyze lists of trajectories (or segments).

    For the most part, this imitates a list. It supports most list behaviors
    (addition, iteration, etc). However, as a list of lists, it contains a
    few useful helpers.

    Parameters
    ----------
    dt : float
        time step between frames in a segment

    Attributes
    ----------
    n_frames : list of int
        number of frames (length) of each segment in this container.
    times : list of float
        duration of each segment (length * self.dt) in this container.
    """
    def __init__(self, segments, dt=None):
        self._segments = segments
        self.dt = dt

    @classmethod
    def from_trajectory_and_indices(cls, trajectory, indices, dt=None):
        segments = [trajectory[idx[0]:idx[1]] for idx in indices]
        return cls(segments, dt)

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


class TrajectoryTransitionAnalysis(object):
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
        self.stateB = transition.stateB - transition.stateA
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
        return {k: self.continuous_segments[k].n_frames
                for k in self.continuous_segments.keys()}

    @property
    def continuous_times(self):
        return {k: self.continuous_segments[k].times
                for k in self.continuous_segments.keys()}

    @property
    def lifetime_frames(self):
        return {k: self.lifetime_segments[k].n_frames
                for k in self.lifetime_segments.keys()}

    @property
    def lifetimes(self):
        return {k: self.lifetime_segments[k].times
                for k in self.lifetime_segments.keys()}

    @property
    def transition_duration_frames(self):
        return {k: self.transition_segments[k].n_frames
                for k in self.transition_segments.keys()}

    @property
    def transition_duration(self):
        return {k: self.transition_segments[k].times
                for k in self.transition_segments.keys()}

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
        return TrajectorySegmentContainer(segments, self.dt)

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
        list of :class:`.Trajectory`
            the frames from (and including) each first entry from `to_vol`
            into `from_vol` until (and including) the next entry into
            `to_vol`, with no frames in `forbidden`, and with frames removed
            from the ends according to `padding`
        """
        if forbidden is None:
            forbidden = paths.EmptyVolume()
        ensemble_BAB = paths.SequentialEnsemble([
            paths.LengthEnsemble(1) & paths.AllInXEnsemble(to_vol),
            paths.AllOutXEnsemble(to_vol) & paths.PartInXEnsemble(from_vol),
            paths.LengthEnsemble(1) & paths.AllInXEnsemble(to_vol)
        ]) & paths.AllOutXEnsemble(forbidden)
        ensemble_AB = paths.SequentialEnsemble([
            paths.LengthEnsemble(1) & paths.AllInXEnsemble(from_vol),
            paths.OptionalEnsemble(paths.AllOutXEnsemble(to_vol)),
            paths.LengthEnsemble(1) & paths.AllInXEnsemble(to_vol)
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

        Returns
        -------
        :class:`.TrajectorySegmentContainer`
            lifetime from first entrance of `stateA` until first entrance of
            `stateB`
        """
        other_state = list(set([self.stateA, self.stateB]) - set([state]))[0]
        segments = self.get_lifetime_segments(
            trajectory=trajectory,
            from_vol=state,
            to_vol=other_state
        )
        return TrajectorySegmentContainer(segments, self.dt)

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

        Returns
        -------
        :class:`.TrajectorySegmentContainer`
            transitions from `stateA` to `stateB` within `trajectory`
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
        return TrajectorySegmentContainer(segments, self.dt)

    def analyze_flux(self, trajectories, state, interface=None):
        """Analysis to obtain flux segments for given state.

        Parameters
        ----------
        trajectories : :class:`.Trajectory` or list of :class:`.Trajectory`
            trajectory to analyze
        state : :class:`.Volume`
            state volume to characterize. Must be one of the states in the
            transition
        interface : :class:`.Volume` or None
            interface to calculate the flux through. If `None`, same as
            `state`

        Returns
        -------
        dict
            key 'in' maps to :class:`.TrajectorySegmentContainer` of frames
            marked "inside" the interface; 'out' maps to frames "outside"
            the interface. The reciprocal of the sum of the mean of these
            two is the flux through the interface.
        """
        if interface is None:
            interface = state
        if isinstance(trajectories, paths.Trajectory):
            trajectories = [trajectories]

        all_flux_dicts = [
            self._analyze_flux_single_traj(traj, state, interface)
            for traj in trajectories
        ]

        empty = TrajectorySegmentContainer([], dt=self.dt)
        total_in = sum([flux['in'] for flux in all_flux_dicts], empty)
        total_out = sum([flux['out'] for flux in all_flux_dicts], empty)
        return {'in': total_in, 'out': total_out}


    def _analyze_flux_single_traj(self, trajectory, state, interface):
        other = list(set([self.stateA, self.stateB]) - set([state]))[0]
        out_segments = self.get_lifetime_segments(
            trajectory=trajectory,
            from_vol=~interface,
            to_vol=state,
            forbidden=other,
            padding=[None, -1]
        )
        out_container = TrajectorySegmentContainer(out_segments, self.dt)
        in_segments = self.get_lifetime_segments(
            trajectory=trajectory,
            from_vol=state,
            to_vol=~interface,
            forbidden=other,
            padding=[None, -1]
        )
        in_container = TrajectorySegmentContainer(in_segments, self.dt)
        return {'in': in_container, 'out': out_container}

    def flux(self, trajectories, state, interface=None):
        if self.dt is None:
            raise RuntimeError("Can't calculate the flux without `dt`")

        flux_dict = self.analyze_flux(trajectories, state, interface)
        return self.flux_from_flux_dict(flux_dict)

    @staticmethod
    def flux_from_flux_dict(flux_dict):
        in_segs = flux_dict['in']
        out_segs = flux_dict['out']
        flux = 1.0 / (np.mean(in_segs.times) + np.mean(out_segs.times))
        return flux

    def analyze(self, trajectories):
        """Full analysis of a trajectory or trajectories.

        Parameters
        ----------
        trajectories : :class:`.Trajectory` or list of :class:`.Trajectory`
        """
        # TODO: I hate using isinstance, but I don't see another way
        if isinstance(trajectories, paths.Trajectory):
            trajectories = [trajectories]

        # shortcuts for readability
        c_segs = self.continuous_segments
        l_segs = self.lifetime_segments
        t_segs = self.transition_segments
        f_dicts = self.flux_segments
        for traj in trajectories:
            for state in [self.stateA, self.stateB]:
                c_segs[state] += self.analyze_continuous_time(traj, state)
                l_segs[state] += self.analyze_lifetime(traj, state)
                f_dict = self.analyze_flux(traj, state)
                f_dicts[state]['in'] += f_dict['in']
                f_dicts[state]['out'] += f_dict['out']
            t_duration_AB = self.analyze_transition_duration(traj,
                                                             self.stateA,
                                                             self.stateB)
            t_duration_BA = self.analyze_transition_duration(traj,
                                                             self.stateB,
                                                             self.stateA)
            t_segs[(self.stateA, self.stateB)] += t_duration_AB
            t_segs[(self.stateB, self.stateA)] += t_duration_BA
        # return self so we can init and analyze in one line
        return self

    # TODO: add a `summary` function to output a nice pandas frame or
    # something
