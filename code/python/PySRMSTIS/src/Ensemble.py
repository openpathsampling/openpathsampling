'''
Created on 03.09.2014

@author: jan-hendrikprinz
'''

class Ensemble(object):
    '''
    
    Example
    -------
    
    >>> EnsembleFactory.TISEnsemble(
    >>>     LambdaVolume(orderparameter_A, 0.0, 0.02), 
    >>>     LambdaVolume(orderparameter_A, 0.0, 0.02), 
    >>>     LambdaVolume(orderparameter_A, 0.0, 0.08), 
    >>>     True
    >>>     )
    '''
    def __init__(self):
        '''
        A path ensemble defines a set of paths
        '''
        
        self._traj = dict()
        
        return
    
    def inside(self, trajectory):
        '''
        Returns True if the trajectory is part of the PathEnsemble
        '''
        return False
    
    def iterate(self, trajectory):
        '''
        Returns true, if the trajectory so far can still be in the ensemble. To check, it assumes that the
        trajectory to length L-1 is okay. This is mainly for interactive usage, when a trajectory is generated.
        
        Notes
        -----
        
        This is only tricky for this that depend on the history like HitXEnsemble or LeaveXEnsembles. In theory these can only
        be checked if the full range of frames has been generated. This could be triggered, when the last frame is reached.
        This is even more difficult if this depends on the length.
        '''

        return False
        
        pass
        
    
    def __str__(self):
        return 'Ensemble'

    def __or__(self, other):
        if self is other:
            return self
        else:
            return EnsembleCombination(self, other, fnc = lambda a,b : a or b, str_fnc = '{0}\nor\n{1}')

    def __xor__(self, other):
        if self is other:
            # A xor A is always Empty
            return EmptyEnsemble()
        else:
            return EnsembleCombination(self, other, fnc = lambda a,b : a ^ b, str_fnc = '{0}\nxor\n{1}')

    def __and__(self, other):
        if self is other:
            return self
        else:
            return EnsembleCombination(self, other, fnc = lambda a,b : a and b, str_fnc = '{0}\nand\n{1}')

    def __sub__(self, other):
        if self is other:
            # A \ A is always Empty
            return EmptyEnsemble()
        else:
            return EnsembleCombination(self, other, fnc = lambda a,b : a and not b, str_fnc = '{0}\nand not\n{1}')
        
    @staticmethod
    def _indent(s):
        spl = s.split('\n')
        spl = ['  ' + p for p in spl]
        return '\n'.join(spl)
    
    def _lencheck(self, trajectory):
        if hasattr(self, 'frames'):
            if type(self.frames) is int:
                return trajectory.frames > self.frames and trajectory.frames >= -self.frames
                
    
class EmptyEnsemble(Ensemble):
    def __init__(self):
        super(EmptyEnsemble, self).__init__()

    def inside(self, trajectory):
        return False

    def iterate(self, trajectory):
        return False

    def __str__(self):
        return 'empty'

class FullEnsemble(Ensemble):
    def __init__(self):
        super(EmptyEnsemble, self).__init__()

    def inside(self, trajectory):
        return True
    
    def iterate(self, trajectory):
        return True

    def __str__(self):
        return 'all'
    
class NegatedEnsemble(Ensemble):
    '''
    Negates an Ensemble and simulates a `not` statement
    '''
    def __init__(self, ensemble):
        super(NegatedEnsemble, self).__init__()
        self.ensemble = ensemble
        
    def inside(self, trajectory):
        return not self.ensemble.inside(trajectory)

    def iterate(self, trajectory):
        return not self.ensemble.iterate(trajectory)

    def __str__(self):
        return 'not ' + str(self.ensemble2)

    
class EnsembleCombination(Ensemble):
    def __init__(self, ensemble1, ensemble2, fnc, str_fnc):
        super(EnsembleCombination, self).__init__()
        self.ensemble1 = ensemble1
        self.ensemble2 = ensemble2
        self.fnc = fnc
        self.sfnc = str_fnc

    def inside(self, trajectory):
        return self.fnc(self.ensemble1.inside(trajectory), self.ensemble2.inside(trajectory))

    def iterate(self, trajectory):
        return self.fnc(self.ensemble1.iterate(trajectory), self.ensemble2.iterate(trajectory))

    def __str__(self):
#        print self.sfnc, self.ensemble1, self.ensemble2, self.sfnc.format('(' + str(self.ensemble1) + ')' , '(' + str(self.ensemble1) + ')')
        return self.sfnc.format('(\n' + Ensemble._indent(str(self.ensemble1)) + '\n)' , '(\n' + Ensemble._indent(str(self.ensemble2)) + '\n)')

class LengthEnsemble(Ensemble):
    def __init__(self, length):
        '''
        A path ensemble that describes path of a specific length
        '''
        
        self.length = length
        pass
    
    def inside(self, trajectory):
        length = trajectory.frames
        if type(self.length) is int:
            return length == self.length
        else:
            return length >= self.length.start and length < self.length.stop
    
    def __str__(self):
        if type(self.length) is int:
            return 'len(x) = {0}'.format(self.length)
        else:
            return 'len(x) in [{0}, {1}]'.format(self.length.start, self.length.stop + 1)

class OutXEnsemble(Ensemble):
    def __init__(self, volume, frames = slice(1,-1), lazy = True):
        '''
        A path ensemble that describes path that start leaves state X in between start and end
        '''
        
        self.volume = volume
        self.frames = frames
        self.lazy = lazy
        pass

    def __str__(self):
        if type(self.frames) is int:
            return 'x[{0}] not in {1}'.format(str(self.frames), str(self.volume))
        else:
            return 'x[t] not in {2} for all t in [{0}:{1}])'.format(self.frames.start, self.frames.stop, self.volume)

    def iterate(self, trajectory):
        pos = trajectory.frames - 1
        checkpoint = -1 
        if type(self.frames) is int:
            if self.frames >= 0:
                if self.frames == pos:
                    checkpoint = pos
        else:
            if self.frames.start >= 0:
                if self.frames.stop >=0:
                    if self.frames.start <= pos and self.frames.stop >= pos:
                        checkpoint = pos
                else:
                    if 1+ pos + self.frames.stop >= self.frames.start:
                        checkpoint = 1 + pos + self.frames.stop
                                        
        if checkpoint >= 0:
#            print 'Out:',checkpoint            
            return not self.volume.inside(trajectory[checkpoint])  
        
        return True    
    
    def inside(self, trajectory):
        if type(self.frames) is int:
            if trajectory.frames > self.frames and trajectory.frames >= -self.frames:
                return not self.volume.inside(trajectory[self.frames])            
            else:
                return False         
        else:
            if self.lazy:
                for s in trajectory[self.frames][0::max(1,trajectory.frames - 1)]:
                    if self.volume.inside(s):
                        return False
            else:
                for s in trajectory[self.frames]:
                    if self.volume.inside(s):
                        return False
                        
            return True

class InXEnsemble(Ensemble):
    def __init__(self, volume, frames, lazy = True):
        '''
        A path ensemble that describes path that start leaves state X in between start and end
        '''
        
        self.volume = volume
        self.frames = frames
        self.lazy = lazy
        pass

    def __str__(self):
        if type(self.frames) is int:
            return 'x[{0}] in {1}'.format(str(self.frames), str(self.volume))
        else:
            return 'x[t] in {2} for all t in [{0}:{1}])'.format(self.frames.start, self.frames.stop, self.volume)

    def iterate(self, trajectory):
        pos = trajectory.frames - 1
        checkpoint = -1 
        if type(self.frames) is int:
            if self.frames >= 0:
                if self.frames == pos:
                    checkpoint = pos
        else:
            if self.frames.start >= 0:
                if self.frames.stop >=0:
                    if self.frames.start <= pos and self.frames.stop >= pos:
                        checkpoint = pos
                else:
                    if 1 + pos + self.frames.stop >= self.frames.start:
                        checkpoint = 1 + pos + self.frames.stop
                                        
        if checkpoint >= 0:
#            print 'In:',checkpoint
            return self.volume.inside(trajectory[checkpoint])  
        
        return True    
    
    
    def inside(self, trajectory):
        '''
        Returns True if the trajectory is part of the PathEnsemble
        '''
                
        if type(self.frames) is int:
            if trajectory.frames > self.frames and trajectory.frames >= -self.frames:
                return self.volume.inside(trajectory[self.frames])            
            else:
                return False
        else:
            if self.lazy:
                for s in trajectory[self.frames][0::max(1,trajectory.frames - 1)]:
                    if not self.volume.inside(s):
                        return False
            else:
                for s in trajectory[self.frames]:
                    if not self.volume.inside(s):
                        return False

            return True
        
        
class HitXEnsemble(Ensemble):
    def __init__(self, volume, frames = slice(1,-1)):
        '''
        A path ensemble that describes path that start leaves state X in between start and end
        '''
        
        self.volume = volume
        self.frames = frames
        pass

    def __str__(self):
        if type(self.frames) is int:
            return 'x[{0}] in {1}'.format(str(self.frames), str(self.volume))
        else:
            return 'x[t] in {2} for one t in [{0}:{1}]'.format(self.frames.start, self.frames.stop, self.volume)


    def iterate(self, trajectory):
        # Always true since a longer trajectory might always make it true
        return True    

    
    def inside(self, trajectory):
        '''
        Returns True if the trajectory is part of the PathEnsemble
        '''
        
        if type(self.frames) is int:
            if trajectory.frames > self.frames and trajectory.frames >= -self.frames:
                return self.volume.inside(trajectory[self.frames])            
            else:
                return False          
        else:
            for s in trajectory[self.frames]:
                if self.volume.inside(s):
                    return True    

            return False

class LeaveXEnsemble(Ensemble):
    def __init__(self, volume, frames = slice(1,-1)):
        '''
        A path ensemble that describes path that start leaves state X in between start and end
        '''
        
        self.volume = volume
        self.frames = frames
        pass
    
    def __str__(self):
        if type(self.frames) is int:
            return 'x[{0}] not in {1}'.format(str(self.frames), str(self.volume))
        else:
            return 'x[t] in {2} for one t in [{0}:{1}])'.format(self.frames.start, self.frames.stop, self.volume)
        
    def iterate(self, trajectory):
        return True
    
    def inside(self, trajectory):
        '''
        Returns True if the trajectory is part of the PathEnsemble
        '''
        
        if type(self.frames) is int:
            if trajectory.frames > self.frames and trajectory.frames >= -self.frames:
                return not self.volume.inside(trajectory[self.frames])            
            else:
                return False       
        else:
            for s in trajectory[self.frames]:
                if not self.volume.inside(s):
                    return True    

            return False

    
class EnsembleFactory():

    @staticmethod
    def StartXEnsemble(volume):
        return InXEnsemble(volume, 0)

    @staticmethod
    def EndXEnsemble(volume):
        return InXEnsemble(volume, -1)

    @staticmethod
    def A2BEnsemble(vol_start, vol_end, lazy = True):
        return (InXEnsemble(vol_start, 0) & InXEnsemble(vol_end, -1)) & OutXEnsemble(vol_start | vol_end, slice(1,-1), lazy)

    @staticmethod
    def TISEnsemble(vol_start, vol_end, vol_interface, lazy = True):
        return (InXEnsemble(vol_start, 0) & InXEnsemble(vol_end, -1)) & (LeaveXEnsemble(vol_interface) & OutXEnsemble(vol_start | vol_end, slice(1,-1), lazy))
    
    hbar = StartXEnsemble



# NOT NEEDED YET
class PathBundle(object):
    def __init__(self):
        '''
        A set of PathEnsembles.
        
        Notes
        -----
        
        Should be handled similar to a set, but with additional was to access different groups of Ensembles
        '''
