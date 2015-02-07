import math
import numpy as np

from openpathsampling.todict import restores_as_full_object
import logging
from ops_logging import initialization_logging
logger = logging.getLogger(__name__)
init_log = logging.getLogger('openpathsampling.initialization')


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

@restores_as_full_object

class ShootingPoint(object):

    def __init__(self, selector, trajectory, index, f = None, sum_bias = None):
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
            The unnormalized probability for the point picked
        sum_bias : float
            The unnormalize probability for the trajectory from which the point was picked

        Notes
        -----
        '''

        self.selector = selector
        self.trajectory = trajectory
        self.index = index
        self._f = f
        self._sum_bias = sum_bias

    @property
    def snapshot(self):
        return self.trajectory[self.index]

    @property
    def sum_bias(self):
        '''
        Return the unnormalized probability for the total trajectory where
        the point has been chosen from.
        
        Notes
        -----
        These partition function like normalizations for a trajectory should
        only be computed only once.  Think about a way to store this. Maybe
        use a cache for the ShootingPoint
        '''
        if self._sum_bias is None:
            self._sum_bias = self.selector.sum_bias(self.trajectory)

        return self._sum_bias
    
    @property
    def probability(self):
        return self._f / self._sum_bias
    
    @property
    def f(self):
        if self._f is None:
            self._f = self.selector.f(self.snapshot, self.trajectory)
            
        return self._f

    @property
    def bias(self):
        return self.f

@restores_as_full_object
class ShootingPointSelector(object):

    @property
    def identifier(self):
        if hasattr(self, 'json'):
            return self.json
        else:
            return None

    def f(self, snapshot, trajectory=None):
        '''
        Returns the unnormalized proposal probability of a snapshot
        
        Notes
        -----
        In principle this is an orderparameter so we could easily add
        caching if useful
        '''
        return 1.0
    
    def probabilities(self, snapshot, trajectory):
        return self.f(snapshot, trajectory) / self.sum_bias(trajectory)
    
    def _biases(self, trajectory):
        '''
        Returns a list of unnormalized proposal probabilities for all
        snapshots in trajectory
        '''
        return [ self.f(s, trajectory) for s in trajectory ]
    
    def sum_bias(self, trajectory):
        '''
        Returns the unnormalized probability probability of a trajectory.
        This is just the sum of all proposal probabilities in a trajectory.
        
        Notes
        -----
        For a uniform distribution this is proportional to the length of the
        trajectory. In this case we can estimate the maximal accepted
        trajectory length for a given acceptance probability.
        
        After we have generated a new trajectory the acceptance probability only for the non-symmetric proposal of
        different snapshots is given by `probability(old_trajectory) / probability(new_trajectory)`
        '''

        return sum(self._biases(trajectory))
    
    def pick(self, trajectory):
        '''
        Returns a ShootingPoint object from which all necessary properties about the selected point can be accessed
        
        Notes
        -----
        
        The native implementation is very slow. Simple picking algorithm should override this function.
        '''
        
        prob_list = self._biases(trajectory)
        sum_bias = sum(prob_list)
        
        rand = np.random.random() * sum_bias
        idx = 0
        prob = prob_list[0]
        while prob <= rand and idx < len(prob_list):
            idx += 1
            prob += prob_list[idx]
            
        point = ShootingPoint(self, trajectory, idx, f = prob_list[idx], sum_bias= sum_bias)

        return point

@restores_as_full_object
class GaussianBiasSelector(ShootingPointSelector):
    def __init__(self, orderparameter, alpha = 1.0, l0 = 0.5):
        '''
        A Selector that biasses according to a specified Orderparameter using a mean l0 and a variance alpha
        '''
        super(GaussianBiasSelector, self).__init__()
        self.orderparameter = orderparameter
        self.alpha = alpha
        self.l0 = l0

    def f(self, snapshot, trajectory=None):
        return math.exp(-self.alpha*(self.orderparameter(snapshot) - self.l0)**2)

@restores_as_full_object
class UniformSelector(ShootingPointSelector):
    def __init__(self, pad_start = 1, pad_end = 1):
        '''
        A Uniform Selector that picks equally likeli any point in the trajectory except the ones excluded using padding
        '''
        super(UniformSelector, self).__init__()
        self.pad_start = pad_start
        self.pad_end = pad_end
        
    def f(self, frame, trajectory=None):
        '''
        Careful, this only returns a correct value for allowed frames since this function does not know about the position in the trajectory
        '''
        return 1.0
    
    def sum_bias(self, trajectory):
        return float(trajectory.frames - self.pad_start - self.pad_end)
        
    def pick(self, trajectory):
        idx = np.random.random_integers(self.pad_start, trajectory.frames - self.pad_end - 1)
        
        point = ShootingPoint(self, trajectory, idx, f = 1.0, sum_bias= self.sum_bias(trajectory))
        
        return point
@restores_as_full_object
class FinalFrameSelector(ShootingPointSelector):
    '''
    Pick final trajectory frame as shooting point.

    This is used for "forward" extension in, e.g., the minus move.
    '''
    def f(self, frame, trajectory):
        if trajectory.index(frame) == len(trajectory)-1:
            return 1.0
        else:
            return 0.0

    def pick(self, trajectory):
        point = ShootingPoint(self, trajectory, len(trajectory)-1, f=1.0, sum_bias=1.0)
        return point

@restores_as_full_object
class FirstFrameSelector(ShootingPointSelector):
    '''
    Pick first trajectory frame as shooting point.

    This is used for "backward" extension in, e.g., the minus move.
    '''
    def f(self, frame, trajectory):
        if trajectory.index(frame) == 0:
            return 1.0
        else:
            return 0.0

    def pick(self, trajectory):
        point = ShootingPoint(self, trajectory, 0, f=1.0, sum_bias=1.0)
        return point
