'''
Created on 07.09.2014

@author: jan-hendrikprinz
'''

class Acceptor(object):
    '''
    An Acceptor is a function that computes an acceptance probability for a given trajectory.
    
    Notes
    -----
    It is only meant to compute acceptance probabilities based on the final proposed trajectory and
    be independent on the proposal probability or the initial trajectory.
    '''


    def __init__(self):
        '''
        Constructor
        '''

    def acceptance(self, trajectory):
        '''
        Returns the acceptance probability for a trajectory
        
        Parameters
        ----------
        
        '''
        return 1.0
    
class InEnsembleAcceptor(Acceptor):
    def __init__(self, ensemble):
        super(InEnsembleAcceptor, self).__init__()
        self.ensemble = ensemble
        
    def acceptance(self, trajectory):
        if self.ensemble.__call__(trajectory):
            return 1.0
        else:
            return 0.0 