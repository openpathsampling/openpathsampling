import math
import numpy as np

class ShootingPointSelector(object):
    def __init__(self):
        pass
    
    def f(self, snapshot):
        '''
        Returns the unnormalized proposal probability of a snapshot
        '''
        return 1.0
    
    def _probability_list(self, trajectory):
        '''
        Returns a list of unnormalized proposal probabilities for all snapshots in trajectory
        '''
        return [ self.f(s) for s in trajectory ]
    
    def bias(self, trajectory):
        '''
        Returns the unnormalized bias probability of a trajectory. This is just the sum of all proposal probabilities in a trajectory.
        
        Notes
        -----
        For a uniform distribution this is proportional to the length of the trajectory.
        In this case we can estimate the maximal accepted trajectory length for a given acceptance probability.
        
        After we have generated a new trajectory the acceptance bias only for the non-symmetric proposal of 
        different snapshots is given by `bias(old_trajectory) / bias(new_trajectory)`
        '''

        return sum(self._probability_list(trajectory))
    
    def pick(self, trajectory):
        '''
        Pick a snapshot from the trajectory according to the picking scheme defined in the actual instance of the ShootingPointSelector
        '''
        idx, bias = self.choose(trajectory)
        return trajectory[idx], bias

    def choose(self, trajectory):
        '''
        Chooses the index of a snapshot from the trajectory according to the picking scheme defined in the actual instance of the ShootingPointSelector
        
        Notes
        -----
        
        The native implementation is very slow. Simple picking algorithm should override this function.
        '''
        
        prob_list = self._probability_list(trajectory)
        bias = sum(prob_list)
        
        rand = np.random.random() * bias
        idx = 0
        prob = prob_list[0]
        while prob <= rand and idx < len(prob_list):
            idx += 1
            prob += prob_list[idx]

        return idx, prob_list[idx] / bias

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
    
    def bias(self, trajectory):
        return trajectory.frames - self.pad_start - self.pad_end
        
    def choose(self, trajectory):
        pos = np.random.random_integers(self.pad_start, trajectory.frames - self.pad_end - 1)
        
        return pos, 1.0 / self.bias(trajectory)