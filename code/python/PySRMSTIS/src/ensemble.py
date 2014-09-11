'''
Created on 03.09.2014

@author: jan-hendrikprinz
'''

class Ensemble(object):
    '''    
    An Ensemble represents a path ensemble, effectively a set of trajectories.
    Typical set operations are allowed, here: and, or, xor, -(without), ~ (inverse = all - x)     
    
    Examples
    --------    
    >>> EnsembleFactory.TISEnsemble(
    >>>     LambdaVolume(orderparameter_A, 0.0, 0.02), 
    >>>     LambdaVolume(orderparameter_A, 0.0, 0.02), 
    >>>     LambdaVolume(orderparameter_A, 0.0, 0.08), 
    >>>     True
    >>>     )

    Notes
    -----
    Maybe replace - by / to get better notation. So war its not been used
    '''
    def __init__(self):
        '''
        A path volume defines a set of paths.
        '''
        
        self._traj = dict()
        self.last = None
        
        return
    
    def __call__(self, trajectory):
        '''
        Returns `True` if the trajectory is part of the path ensemble.
        '''
        return False
    
    def forward(self, trajectory):
        '''
        Returns true, if the trajectory so far can still be in the ensemble if it is appended by a frame. To check, it assumes that the
        trajectory to length L-1 is okay. This is mainly for interactive usage, when a trajectory is generated.
        
        Parameters
        ----------
        trajectory : Trajectory
            the actual trajectory to be tested
        
        Returns
        -------
        forward : bool
            Returns true or false if using a forward step (extending the trajectory forward in time at its end) `trajectory` could  
            still be in the ensemble and thus makes sense to continue a simulation
        

        Notes
        -----
        This is only tricky for this that depend on the history like HitXEnsemble or LeaveXEnsembles. In theory these can only
        be checked if the full range of frames has been generated. This could be triggered, when the last frame is reached.
        This is even more difficult if this depends on the length.
        '''

        return False        
        pass
    
    def backward(self, trajectory):
        '''
        Returns true, if the trajectory so far can still be in the ensemble if it is prepended by a frame. To check, it assumes that the
        trajectory from index 1 is okay. This is mainly for interactive usage, when a trajectory is generated using a backward move.
        
        Parameters
        ----------
        trajectory : Trajectory
            the actual trajectory to be tested
        
        Returns
        -------
        backward : bool
            Returns true or false if using a backward step (extending the trajectory backwards in time at its beginning) `trajectory` could  
            still be in the ensemble and thus makes sense to continue a simulation
        
        Notes
        
        This is only tricky for this that depend on the history like HitXEnsemble or LeaveXEnsembles. In theory these can only
        be checked if the full range of frames has been generated. This could be triggered, when the last frame is reached.
        This is even more difficult if this depends on the length.
        '''

        return False        
        pass

    def __str__(self):
        '''
        Returns a complete mathematical expression that defines the current ensemble in a readable form.
        
        Notes
        -----
        This should be cleaned up a little
        '''
        return 'Ensemble'

    def __or__(self, other):
        if self is other:
            return self
        elif type(self) is EmptyEnsemble:
            return other
        elif type(other) is EmptyEnsemble:
            return self
        elif type(self) is FullEnsemble:
            return self
        elif type(other) is FullEnsemble:
            return other
        else:
            return EnsembleCombination(self, other, fnc = lambda a,b : a or b, str_fnc = '{0}\nor\n{1}')

    def __xor__(self, other):
        if self is other:
            return EmptyEnsemble()
        elif type(self) is EmptyEnsemble:
            return other
        elif type(other) is EmptyEnsemble:
            return self
        elif type(self) is FullEnsemble:
            return NegatedEnsemble(other)
        elif type(other) is FullEnsemble:
            return NegatedEnsemble(self)        
        else:
            return EnsembleCombination(self, other, fnc = lambda a,b : a ^ b, str_fnc = '{0}\nxor\n{1}')

    def __and__(self, other):
        if self is other:
            return self
        elif type(self) is EmptyEnsemble:
            return self
        elif type(other) is EmptyEnsemble:
            return other
        elif type(self) is FullEnsemble:
            return other
        elif type(other) is FullEnsemble:
            return self
        else:
            return EnsembleCombination(self, other, fnc = lambda a,b : a and b, str_fnc = '{0}\nand\n{1}')

    def __sub__(self, other):
        if self is other:
            return EmptyEnsemble()
        elif type(self) is EmptyEnsemble:
            return self
        elif type(other) is EmptyEnsemble:
            return self
        elif type(self) is FullEnsemble:
            return NegatedEnsemble(other)
        elif type(other) is FullEnsemble:
            return EmptyEnsemble()        
        else:
            return EnsembleCombination(self, other, fnc = lambda a,b : a and not b, str_fnc = '{0}\nand not\n{1}')
        
    def __invert__(self):
        if type(self) is EmptyEnsemble:
            return FullEnsemble()
        elif type(self) is FullEnsemble:
            return EmptyEnsemble(self)
        else:
            return NegatedEnsemble(self)

        
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
    '''
    The empty path ensemble of no trajectories.
    '''
    def __init__(self):
        super(EmptyEnsemble, self).__init__()

    def __call__(self, trajectory):
        return False

    def forward(self, trajectory):
        return False

    def backward(self, trajectory):
        return False

    def __str__(self):
        return 'empty'

class FullEnsemble(Ensemble):
    '''
    The full path ensemble of all possible trajectories.
    '''
    def __init__(self):
        super(EmptyEnsemble, self).__init__()

    def __call__(self, trajectory):
        return True
    
    def forward(self, trajectory):
        return True

    def backward(self, trajectory):
        return True

    def __str__(self):
        return 'all'
    
class NegatedEnsemble(Ensemble):
    '''
    Negates an Ensemble and simulates a `not` statement
    '''
    def __init__(self, volume):
        super(NegatedEnsemble, self).__init__()
        self.ensemble = volume
        
    def __call__(self, trajectory):
        return not self.ensemble(trajectory)

    def forward(self, trajectory):
        return not self.ensemble.forward(trajectory)

    def backward(self, trajectory):
        return not self.ensemble.forward(trajectory)

    def __str__(self):
        return 'not ' + str(self.ensemble2)

    
class EnsembleCombination(Ensemble):
    '''
    Represent the boolean concatenation of two ensembles
    '''
    def __init__(self, ensemble1, ensemble2, fnc, str_fnc):
        super(EnsembleCombination, self).__init__()
        self.ensemble1 = ensemble1
        self.ensemble2 = ensemble2
        self.fnc = fnc
        self.sfnc = str_fnc

    def __call__(self, trajectory):
        return self.fnc(self.ensemble1(trajectory), self.ensemble2(trajectory))

    def forward(self, trajectory):
        return self.fnc(self.ensemble1.forward(trajectory), self.ensemble2.forward(trajectory))

    def backward(self, trajectory):
        return self.fnc(self.ensemble1.backward(trajectory), self.ensemble2.backward(trajectory))

    def __str__(self):
#        print self.sfnc, self.ensemble1, self.ensemble2, self.sfnc.format('(' + str(self.ensemble1) + ')' , '(' + str(self.ensemble1) + ')')
        return self.sfnc.format('(\n' + Ensemble._indent(str(self.ensemble1)) + '\n)' , '(\n' + Ensemble._indent(str(self.ensemble2)) + '\n)')

class LengthEnsemble(Ensemble):
    '''
    Represents an ensemble the contains trajectories of a specific length
    '''
    def __init__(self, length):
        '''
        A path ensemble that describes path of a specific length
        
        Parameters
        ----------
        length : int or slice
            The specific length (int) or the range of allowed trajectory lengths (slice)
        '''
        
        self.length = length
        pass
    
    def __call__(self, trajectory):
        length = trajectory.frames
        if type(self.length) is int:
            return length == self.length
        else:
            return length >= self.length.start and length < self.length.stop
        
    def forward(self, trajectory):
        length = trajectory.frames
        if type(self.length) is int:
            return length < self.length
        else:
            return length < self.length.stop

    def backward(self, trajectory):
        length = trajectory.frames
        if type(self.length) is int:
            return length < self.length
        else:
            return length < self.length.stop
                
    
    def __str__(self):
        if type(self.length) is int:
            return 'len(x) = {0}'.format(self.length)
        else:
            return 'len(x) in [{0}, {1}]'.format(self.length.start, self.length.stop + 1)
        
class VolumeEnsemble(Ensemble):
    '''
    Describes an path ensemble using a volume object
    '''    
    def __init__(self, volume, frames = slice(1,-1), lazy = True):
        super(VolumeEnsemble, self).__init__()
        self._volume = volume
        self.frames = frames
        self.lazy = lazy
        pass
    
    @property
    def volume(self):
        '''
        The volume that is used in the specification.
        '''
        return self._volume
    
class InXEnsemble(VolumeEnsemble):
    '''
    Represents an ensemble where part of the trajectory is in a specified volume
    '''
    
    def forward(self, trajectory):
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
#            print 'Out:',checkpoint          
            return self.volume(trajectory[checkpoint])  
        
        return True    
    
    def backward(self, trajectory):
        pos = - trajectory.frames
        checkpoint = -1
        if type(self.frames) is int:
            if self.frames < 0:
                if self.frames == pos:
                    checkpoint = pos
        else:
            if self.frames.stop < 0:
                if self.frames.start < 0:
                    if self.frames.start <= pos and self.frames.stop >= pos:
                        checkpoint = pos
                else:
                    if pos + self.frames.start - 1 <= self.frames.stop:
                        checkpoint = self.frames.start - 1
                                        
        if checkpoint >= 0:
#            print 'Out:',checkpoint            
            return self.volume(trajectory[checkpoint])  
        
        return True    
    
    def __call__(self, trajectory):
        if type(self.frames) is int:
            if trajectory.frames > self.frames and trajectory.frames >= -self.frames:
                return self.volume(trajectory[self.frames])            
            else:
                return True
        else:
            if self.lazy:
                for s in trajectory[self.frames][0::max(1,trajectory.frames - 1)]:
                    if not self.volume(s):
                        return False
            else:
                for s in trajectory[self.frames]:
                    if not self.volume(s):
                        return False
                        
            return True


class OutXEnsemble(InXEnsemble):
    '''
    Represents an ensemble where part of the trajectory is outside a specified volume
    '''    
    @property
    def volume(self):
        return ~ self._volume
    
    def __str__(self):
        if type(self.frames) is int:
            return 'x[{0}] not in {1}'.format(str(self.frames), str(self.volume))
        else:
            return 'x[t] not in {2} for all t in [{0}:{1}])'.format(self.frames.start, self.frames.stop, self.volume)



        
class HitXEnsemble(VolumeEnsemble):
    '''
    Represents an ensemble where part of the trajectory visits a specified volume at least for one frame
    '''

    def __str__(self):
        if type(self.frames) is int:
            return 'x[{0}] in {1}'.format(str(self.frames), str(self.volume))
        else:
            return 'x[t] in {2} for one t in [{0}:{1}]'.format(self.frames.start, self.frames.stop, self.volume)


    def forward(self, trajectory):
        pos = trajectory.frames - 1
        if type(self.frames) is int:
            if self.frames >= 0:
                if self.frames == pos:
                    return self.volume(trajectory[pos])  
        else:
            if self.lazy:
                # There is no way to test in a lazy way without storing the result from previous frames. So keep on running
                return True
            else:
                if self.frames.stop >= 0:
                    if self.frames.start >= 0:
                        if 1 + self.frames.stop <= pos:
                            return self(trajectory)
                    else:
                        return True
                                        
        return True    

    def backward(self, trajectory):
        pos = - trajectory.frames
        if type(self.frames) is int:
            if self.frames < 0:
                if self.frames == pos:
                    return self.volume(trajectory[pos])  
        else:
            if self.lazy:
                # There is no way to test in a lazy way without storing the result from previous frames
                return True
            else:
                if self.frames.stop < 0:
                    if self.frames.start < 0:
                        if self.frames.start >= pos:
                            return self(trajectory)
                    else:
                        return True
                                                
        return True    

    
    def __call__(self, trajectory):
        '''
        Returns True if the trajectory is part of the PathEnsemble
        
        Parameters
        ----------
        trajectory : Trajectory
            The trajectory to be checked
        '''
        
        if type(self.frames) is int:
            if trajectory.frames > self.frames and trajectory.frames >= -self.frames:
                return self.volume(trajectory[self.frames])            
            else:
                return False          
        else:
            for s in trajectory[self.frames]:
                if self.volume(s):
                    return True    

            return False

class LeaveXEnsemble(HitXEnsemble):
    '''
    Represents an ensemble where part of the trajectory leaves a specified volume at least for one frame
    '''
    def __str__(self):
        if type(self.frames) is int:
            return 'x[{0}] not in {1}'.format(str(self.frames), str(self.volume))
        else:
            return 'x[t] in {2} for one t in [{0}:{1}])'.format(self.frames.start, self.frames.stop, self.volume)
        
    @property
    def volume(self):
        return ~ self._volume
        
class AlteredTrajectoryEnsemble(Ensemble):
    '''
    Represents an ensemble where an altered version of a trajectory (extended, reversed, cropped) is part of a given ensemble
    '''
    def __init__(self, ensemble):
        '''
        Represents an ensemble which is the given ensemble but for trajectories where some trajectory is prepended
        '''
        
        super(AlteredTrajectoryEnsemble, self).__init__()
        self.ensemble = ensemble
                
    def _alter(self, trajectory):
        return trajectory
        
    def __call__(self, trajectory):
        return self.ensemble(self, self._alter(trajectory))

    def forward(self, trajectory):
        return self.ensemble.forward(self._alter(trajectory))

    def backward(self, trajectory):
        return self.ensemble.backward(self._alter(trajectory))

class BackwardPrependedTrajectoryEnsemble(AlteredTrajectoryEnsemble):
    '''
    Represents an ensemble which is the given ensemble but for trajectories where some trajectory is prepended
    '''
    def __init__(self, ensemble, trajectory):        
        super(BackwardPrependedTrajectoryEnsemble, self).__init__(ensemble)
        self.add_traj = trajectory        

    def _alter(self, trajectory):
#        print [ s.idx for s in trajectory.reversed + self.add_traj]
        return trajectory.reversed + self.add_traj

class ForwardAppendedTrajectoryEnsemble(AlteredTrajectoryEnsemble):
    '''
    Represents an ensemble which is the given ensemble but for trajectories where some trajectory is appended
    '''
    def __init__(self, ensemble, trajectory):
        super(ForwardAppendedTrajectoryEnsemble, self).__init__(ensemble)
        self.add_traj = trajectory

    def _alter(self, trajectory):
        return self.add_traj + trajectory
    
class ReversedTrajectoryEnsemble(AlteredTrajectoryEnsemble):
    '''
    Represents an ensemble 
    '''
    def _alter(self, trajectory):
        return trajectory.reverse()
    
class EnsembleFactory():
    '''
    Convenience class to construct Ensembles
    '''
    @staticmethod
    def StartXEnsemble(volume):
        '''
        Construct an ensemble that starts (x[0]) in the specified volume
        
        Parameters
        ----------
        volume : volume
            The volume to start in 
        
        Returns
        -------
        ensemble : Ensemble
            The constructed Ensemble
        '''
        return InXEnsemble(volume, 0)

    @staticmethod
    def EndXEnsemble(volume):
        '''
        Construct an ensemble that ends (x[-1]) in the specified volume
        
        Parameters
        ----------
        volume : volume
            The volume to end in 
        
        Returns
        -------
        ensemble : Ensemble
            The constructed Ensemble
        '''        
        return InXEnsemble(volume, -1)

    @staticmethod
    def A2BEnsemble(volume_a, volume_b, lazy = True):
        '''
        Construct an ensemble that starts in (x[0]) in volume_a, ends in volume_b and is in either volumes in between
        
        Parameters
        ----------
        volume_a : volume
            The volume to start in 
        volume_b : volume
            The volume to end in 
        
        Returns
        -------
        ensemble : Ensemble
            The constructed Ensemble
        '''        
        return (InXEnsemble(volume_a, 0) & InXEnsemble(volume_b, -1)) & OutXEnsemble(volume_a | volume_b, slice(1,-1), lazy)

    @staticmethod
    def TISEnsemble(volume_a, volume_b, volume_x, lazy = True):
        '''
        Construct an TIS ensemble that starts in (x[0]) in volume_a, ends in volume_b and is in either volumes in between
        and will also leave volume_x at some point
        
        Parameters
        ----------
        volume_a : volume
            The volume to start in 
        volume_b : volume
            The volume to end in 
        volume_x : volume
            The volume to leave 
        
        Returns
        -------
        ensemble : Ensemble
            The constructed Ensemble
        '''        

        return (InXEnsemble(volume_a, 0) & InXEnsemble(volume_b, -1)) & (LeaveXEnsemble(volume_x) & OutXEnsemble(volume_a | volume_b, slice(1,-1), lazy))
    
    hbar = StartXEnsemble