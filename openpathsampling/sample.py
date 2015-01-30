import random

import openpathsampling as paths
import copy

class SampleKeyError(Exception):
    def __init__(self, key, sample, sample_key):
        self.key = key
        self.sample = sample
        self.sample_key = sample_key
        self.msg = (str(self.key) + " does not match " + str(self.sample_key)
                    + " from " + str(self.sample))

class Sample(object):
    """
    A Sample represents a given "draw" from its ensemble, and is the return
    object from a PathMover. It and contains all information about the move,
    initial trajectories, new trajectories (both as references). 
    
    Since each Sample is a single representative of a single ensemble, each
    Sample consists of one replica ID, one trajectory, and one ensemble.
    This means that movers which generate more than one "draw" (often from
    different ensembles, e.g. replica exchange) will generate more than one
    Sample object.

    Attributes
    ----------
    replica : integer
        The replica ID to which this Sample applies
    trajectory : Trajectory
        The trajectory (path) for this sample
    ensemble : Ensemble
        The Ensemble this sample is drawn from
    details : MoveDetails
        Object 
    step : integer
        the Monte Carlo step number associated with this Sample
    """

    def __init__(self, replica=None, trajectory=None, ensemble=None, details=None, step=-1):
        self.replica = replica
        self.ensemble = ensemble
        self.trajectory = trajectory
        self.details = details
        self.step = step

    def __call__(self):
        return self.trajectory

    def __str__(self):
        mystr = "Step: "+str(self.step)+"\n"
        mystr += "Replica: "+str(self.replica)+"\n"
        mystr += "Trajectory: "+str(self.trajectory)+"\n"
        mystr += "Ensemble: "+repr(self.ensemble)+"\n"
        mystr += "Details: "+str(self.details)+"\n"
        return mystr

    @staticmethod
    def set_time(step, samples):
        for sample in samples:
            sample.step = step
