import openpathsampling as paths
import numpy as np

# TODO: I think this should be a method in Trajectory
def subtrajectory_indices(trajectory, subtrajectories):
    results = []
    for subtraj in subtrajectories:
        results.append([trajectory.index(s) for s in subtraj])
    return results


class SingleTrajectoryAnalysis(object):
    def __init__(self, transition, dt=None):
        self.transition = transition
        self.dt = dt
        self.stateA = transition.stateA
        self.stateB = transition.stateB
        self.continuous_time = {self.stateA : [], self.stateB : []}
        self.lifetime = {self.stateA : [], self.stateB : []}
        self.flux = {self.stateA : {}, self.stateB : {}}


    def analyze_continuous_time(self, trajectory, state):
        ensemble = paths.AllInXEnsemble(state)
        segments = ensemble.split(trajectory)
        lengths = [len(seg) for seg in segments]
        self.continuous_time[state] += lengths

    def analyze_lifetime(self, trajectory, state):
        other_state = list(set([self.stateA, self.stateB]) - set([state]))[0]
	ensemble_BAB = paths.SequentialEnsemble([
	    paths.AllInXEnsemble(other_state) & paths.LengthEnsemble(1),
	    paths.PartInXEnsemble(state) & paths.AllOutXEnsemble(other_state),
	    paths.AllInXEnsemble(other_state) & paths.LengthEnsemble(1)
	])
	ensemble_AB = paths.SequentialEnsemble([
	    paths.AllInXEnsemble(state) & paths.LengthEnsemble(1),
	    paths.OptionalEnsemble(
		paths.PartInXEnsemble(state) &
		paths.AllOutXEnsemble(other_state)
	    ),
	    paths.AllInXEnsemble(other_state) & paths.LengthEnsemble(1)
	])


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

    def continuous_time_distribution(self, state):
        dist = np.array(self.continuous_time[state])
        if self.dt is not None:
            dist *= self.dt
        return dist

    def lifetime_distribution(self, state):
        dist = np.array(self.lifetime[state])
        if self.dt is not None:
            dist *= self.dt
        return dist

    def summary(self):
        # return a nice text/pandas output
        pass
