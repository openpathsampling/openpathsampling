'''
Created on 06.09.2014

@author: jan-hendrikprinz
'''

class Stopper(object):
    '''
    A Stopper works on trajectories and/or snapshots and is called by the generator and
    returns True if the generator should be stop and no more frames should be generated. Usually these are generated from Ensembles
    or Volumes.
    '''


    def __init__(self, params):
        '''
        Constructor
        '''
        
        self.trajectory = None
        self.status = None
        
        pass
        
    def __call__(self, trajectory):
        '''
        If not properly defined always stop the simulation
        '''

        self.trajectory = trajectory
        self.status = self._stop(trajectory)
        return self.status
    
    
class MaxLengthStopper(Stopper):
    '''
    A Stopper that stops a generator when the trajectory has reached a maximal length
    '''
    def __init__(self, length):
        self.length = length
        
    def _stop(self, trajectory):
        return len(trajectory) < self.length
    
class EnsembleStopper(Stopper):
    '''
    A Stopper that stops when the generated trajectory cannot be in the specified ensemble anymore.
    This is based on Ensemble.forward. Make sure that this is really possible to achieve. E.g.
    A trajectory that should end in a specific state can be arbitrarily long so that the stopper is
    never triggered.
    '''
    def __init__(self, ensemble):
        self.ensemble = ensemble
        
    def _stop(self, trajectory):
        return self.ensemble.forward(trajectory)

class VolumeStopper(Stopper):
    '''
    A Stopper that stops when the specified Volume is hit. trajectory[-1] in Volume
    '''
    def __init__(self, volume):
        self.volume = volume
        
    def _stop(self, trajectory):
        return self.volume.outside(trajectory[-1])