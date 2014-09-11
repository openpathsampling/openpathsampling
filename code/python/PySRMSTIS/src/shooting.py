import math
import numpy as np

#############################################################################
#
#
#
#
# Notes
# -----
# 
#
#  
#############################################################################

class ShootingPoint(object):
    def __init__(self, selector, trajectory, index, f = None, Z = None):
        '''
        Constructs a ShootingPoint object.
        
        parameters
        ----------
        
        selector : ShootingPointSelector()
            The Selector used to generate the seleted point
        trajectory : Trajectory()
            The parent trajectory a point is selected from.
        index : int
            The actual index of the point picked from the trajectory
        f : float
            The unnormalized bias for the point picked
        Z : float 
            The unnormalize bias for the trajectory from which the point was picked

        Notes
        -----
        '''

        self.selector = selector
        self.trajectory = trajectory
        self.idx = index
        self._f = None
        self._Z = None
        
    @property
    def snapshot(self):
        return self.trajectory[self.idx]

    @property
    def Z(self):
        '''
        Return the unnormalized bias for the total trajectory where the point has been chosen from.
        
        Notes
        -----
        These partition function like normalizations for a trajectoy should only be computed only once.
        Think about a way to store this. Maybe use a cache for the ShootingPoint
        '''
        if self._Z is None:
            self._Z = self.selector.Z(self.trajectory)

        return self._Z
    
    @property
    def bias(self):
        return self._f / self._Z
    
    @property
    def f(self):
        if self._f is None:
            self._f = self.selector.f(self.snapshot)
            
        return self._f
    
    @property
    def index(self):
        return self.idx
    
class ShootingPointSelector(object):
    def __init__(self):
        pass
    
    def f(self, snapshot):
        '''
        Returns the unnormalized proposal probability of a snapshot
        
        Notes
        -----
        
        In principle this is an orderparameter so we could easily add caching if useful
        '''
        return 1.0
    
    def bias(self, snapshot, trajectory):
        return self.f(snapshot) / self.Z(trajectory)
    
    def _probabilities(self, trajectory):
        '''
        Returns a list of unnormalized proposal probabilities for all snapshots in trajectory
        '''
        return [ self.f(s) for s in trajectory ]
    
    def Z(self, trajectory):
        '''
        Returns the unnormalized bias probability of a trajectory. This is just the sum of all proposal probabilities in a trajectory.
        
        Notes
        -----
        For a uniform distribution this is proportional to the length of the trajectory.
        In this case we can estimate the maximal accepted trajectory length for a given acceptance probability.
        
        After we have generated a new trajectory the acceptance bias only for the non-symmetric proposal of 
        different snapshots is given by `bias(old_trajectory) / bias(new_trajectory)`
        '''

        return sum(self._probabilities(trajectory))
    
    def pick(self, trajectory):
        '''
        Returns a ShootingPoint object from which all necessary properties about the selected point can be accessed
        
        Notes
        -----
        
        The native implementation is very slow. Simple picking algorithm should override this function.
        '''
        
        prob_list = self._probabilities(trajectory)
        Z = sum(prob_list)
        
        rand = np.random.random() * Z
        idx = 0
        prob = prob_list[0]
        while prob <= rand and idx < len(prob_list):
            idx += 1
            prob += prob_list[idx]
            
        point = ShootingPoint(self, trajectory, idx, f = prob_list[idx], Z = Z)

        return point

class GaussianBiasSelector(ShootingPointSelector):
    def __init__(self, orderparameter, alpha = 1.0, l0 = 0.5):
        '''
        A Selector that biasses according to a specified Orderparameter using a mean l0 and a variance alpha
        '''
        super(GaussianBiasSelector, self).__init__()
        self.orderparameter = orderparameter
        self.alpha = alpha
        self.l0 = l0

    def f(self, snapshot):
        return math.exp(-self.alpha*(self.orderparameter(snapshot) - self.l0)**2)

class UniformSelector(ShootingPointSelector):
    def __init__(self, pad_start = 1, pad_end = 1):
        '''
        A Uniform Selector that picks equally likeli any point in the trajectory except the ones excluded using padding
        '''
        super(UniformSelector, self).__init__()
        self.pad_start = pad_start
        self.pad_end = pad_end
        
    def f(self, frame):
        '''
        Careful, this only returns a correct value for allowed frames since this function does not know about the position in the trajectory
        '''
        return 1.0
    
    def Z(self, trajectory):
        return float(trajectory.frames - self.pad_start - self.pad_end)
        
    def pick(self, trajectory):
        idx = np.random.random_integers(self.pad_start, trajectory.frames - self.pad_end - 1)
        
        point = ShootingPoint(self, trajectory, idx, f = 1.0, Z = self.Z(trajectory))
        
        return point