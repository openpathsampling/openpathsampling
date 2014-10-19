'''
Created on 03.09.2014

@author: jan-hendrikprinz, David W.H. Swenson
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

    use_shortcircuit = True

    def __init__(self):
        '''
        A path volume defines a set of paths.
        '''
        
        self._traj = dict()
        self.last = None
        
        return
    
    def __call__(self, trajectory, lazy=None):
        '''
        Returns `True` if the trajectory is part of the path ensemble.

        Parameters
        ----------
        lazy : boolean
            If lazy is not None it overrides the default setting in the ensemble
        '''
        return False

    def check(self, trajectory):
        return self(trajectory, lazy = False)
    
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

        return True        
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

        return True        
        pass

    can_append = forward
    can_prepend = backward


    def locate(self, trajectory, lazy=True, max_length=None, min_length=1, overlap=1):
        '''
        Returns a list of trajectories that contain sub-trajectories which are in the given ensemble.

        Parameters
        ----------
        trajectory : Trajectory
            the actual trajectory to be splitted into ensemble parts
        lazy : boolean
            if True will use a faster almost linear algorithm, while False will run through all possibilities starting with the
            largest ones
        max_length : int > 0
            if set this determines the maximal size to be tested (is mainly used in the recursion)
        min_length : int > 0
            if set this determines the minimal size to be tested (in lazy mode might no
        overlap : int >= 0
            determines the allowed overlap of all trajectories to be found. A value of x means that
            two sub-trajectorie can share up to x frames at the beginning and x frames at the end.
            Default is 1

        Returns
        -------
        list of slices
            Returns a list of index-slices for sub-trajectories in trajectory that are in the ensemble.
        '''
        ensemble_list = []

        length = len(trajectory)

        if max_length is None:
            max_length = length

        max_length = min(length, max_length)
        min_length = max(1, min_length)

        if not lazy:
            # this tries all possible sub-trajectories starting with the longest ones and uses recursion
            for l in range(max_length,min_length - 1,-1):
                for start in range(0,length-l+1):
                    tt = trajectory[start:start+l]
#                    print start, start+l
                    # test using lazy=False
                    if self(tt, lazy=False):
                        list_left = []
                        list_right = []
                        if l > min_length:
                            pad = min(overlap, l - 1)
                            tt_left = trajectory[0:start + pad]
                            list_left = self.locate(tt_left, max_length=l)

                            tt_right = trajectory[start + l - pad:length]
                            list_right = self.locate(tt_right, max_length=l)

#                        ensemble_list = list_left + [tt] + list_right
                        ensemble_list = list_left + [slice(start,start+l)] + list_right

                        # no need to look further inside the iterations caught everything!
                        break
                else:
                    continue
                break

            return ensemble_list
        else:
            start = 0
            end = min_length

            while start <= length - min_length and end <= length:
                tt = trajectory[start:end]
#                print start,end
                if self.forward(tt) and end<length:
                    end += 1
                else:
                    if self(tt, lazy=False):
                        ensemble_list.append(slice(start,end))
                        pad = min(overlap, end - start - 1)
                        start = end - pad
                    else:
                        start += 1
                    end = start + min_length

            return ensemble_list

    def split(self, trajectory, lazy=True, max_length=None, min_length=1, overlap=1):
        '''
        Returns a list of trajectories that contain sub-trajectories which are in the given ensemble.

        Parameters
        ----------
        trajectory : Trajectory
            the actual trajectory to be splitted into ensemble parts
        lazy : boolean
            if True will use a faster almost linear algorithm, while False will run through all possibilities starting with the
            largest ones
        max_length : int > 0
            if set this determines the maximal size to be tested (is mainly used in the recursion)
        min_length : int > 0
            if set this determines the minimal size to be tested (in lazy mode might no
        overlap : int >= 0
            determines the allowed overlap of all trajectories to be found. A value of x means that
            two sub-trajectory can share up to x frames at the beginning and x frames at the end.
            Default is 1

        Returns
        -------
        list of Trajectory
            Returns a list of sub-trajectories in trajectory that are in the ensemble.

        Notes
        -----
        This uses self.locate and returns the actual sub-trajectories
        '''

        indices = self.locate(trajectory, lazy, max_length, min_length, overlap)

        return [trajectory[part] for part in indices]


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

    def __call__(self, trajectory, lazy=None):
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
        super(FullEnsemble, self).__init__()

    def __call__(self, trajectory, lazy=None):
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
        
    def __call__(self, trajectory, lazy=None):
        return not self.ensemble(trajectory, lazy)

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

    def __call__(self, trajectory, lazy=None):
        # Shortcircuit will automatically skip the second part of the combination if the result does not depend on it!
        # This makes sense since the expensive part is the ensemble testing not computing two logic operations
        if Ensemble.use_shortcircuit:
            a = self.ensemble1(trajectory, lazy)
            res_true = self.fnc(a, True)
            res_false = self.fnc(a, False)
            if res_false == res_true:
                # result is independent of ensemble_b so ignore it
                return res_true
            else:
                b = self.ensemble2(trajectory, lazy)
                return self.fnc(a, b)
        else:
            return self.fnc(self.ensemble1(trajectory, lazy), self.ensemble2(trajectory, lazy))

    def forward(self, trajectory):
        if Ensemble.use_shortcircuit:
            a = self.ensemble1.can_append(trajectory)
            res_true = self.fnc(a, True)
            res_false = self.fnc(a, False)
            if res_false == res_true:
                # result is independent of ensemble_b so ignore it
                return res_true
            else:
                b = self.ensemble2.can_append(trajectory)
                return self.fnc(a, b)
        else:
            return self.fnc(self.ensemble1.can_append(trajectory), self.ensemble2.can_append(trajectory))

    def backward(self, trajectory):
        if Ensemble.use_shortcircuit:
            a = self.ensemble1.backward(trajectory)
            res_true = self.fnc(a, True)
            res_false = self.fnc(a, False)
            if res_false == res_true:
                # result is independent of ensemble_b so ignore it
                return res_true
            else:
                b = self.ensemble2.backward(trajectory)
                return self.fnc(a, b)
        else:
            return self.fnc(self.ensemble1.backward(trajectory), self.ensemble2.backward(trajectory))

    can_append = forward
    can_prepend = backward

    def __str__(self):
#        print self.sfnc, self.ensemble1, self.ensemble2, self.sfnc.format('(' + str(self.ensemble1) + ')' , '(' + str(self.ensemble1) + ')')
        return self.sfnc.format('(\n' + Ensemble._indent(str(self.ensemble1)) + '\n)' , '(\n' + Ensemble._indent(str(self.ensemble2)) + '\n)')

class SequentialEnsemble(Ensemble):
    """
    An ensemble that consists of several ensembles which are satisfied by
    the trajectory in sequence. 

    Attributes
    ----------
    ensembles : tuple of Ensemble
        The ensembles, in time-order of when they should occur in the
        trajectory.
    min_overlap : int or tuple of int
        The minimum number of frames that overlap between two ensembles in
        the sequence. A positive number n indicates that at least n frames
        must be in both ensembles at the transition between them. A negative
        number -n indicates that at least n frames in neither ensemble at
        the transition between them. If given as a list, the list should be
        of length len(ensembles)-1, with one value for each transition. If
        given as an integer, that value will be used for all transitions.
    max_overlap : int or list of int
        The maximum number of frames that overlap between two ensembles in
        the sequence. A positive number n indicates that no more than n
        frames can be in both ensembles at the transition between them. A
        negative number -n indicates no more than n frames in neither
        ensemble at the transition between them. If given as a list, the
        list should be of length len(ensembles)-1, with one value for each
        transition. If given as an integer, that value will be used for all
        transitions.

    Todo
    ----
        Overlap features not implemented because ohmygod this was hard
        enough already.
    """

    def __init__(self, ensembles, min_overlap=0, max_overlap=0, greedy=False):
        # make tuples of the min/max overlaps
        if type(min_overlap) is int:
            min_overlap = (min_overlap, )*(len(ensembles)-1)
        if type(max_overlap) is int:
            max_overlap = (max_overlap, )*(len(ensembles)-1)

        self.ensembles = tuple(ensembles)
        self.min_overlap = tuple(min_overlap)
        self.max_overlap = tuple(max_overlap)
        self.greedy = greedy

        # sanity checks
        if len(self.min_overlap) != len(self.max_overlap):
            raise ValueError("len(min_overlap) != len(max_overlap)")
        if len(self.min_overlap) != len(self.ensembles)-1:
            raise ValueError("Number of overlaps doesn't match number of transitions")
        for i in range(len(self.min_overlap)):
            if min_overlap[i] > max_overlap[i]:
                raise ValueError("min_overlap greater than max_overlap!")


    def __call__(self, trajectory, lazy=None):
        # it is easiest to understand this decision tree as a simplified
        # version of the can_append decision tree; see that for detailed
        # comments
        ens_num = 0
        subtraj_first = 0
        traj_final = len(trajectory)
        final_ens = len(self.ensembles)-1
        while True:
            subtraj_final = self._find_subtraj_final(trajectory, 
                                                     subtraj_first, ens_num)
            if subtraj_final - subtraj_first > 0:
                subtraj = trajectory[slice(subtraj_first, subtraj_final)]
                if ens_num == final_ens:
                    if subtraj_final == traj_final:
                        return True # assigned all frames; all ensembles
                    else:
                        return False
                else:
                    ens_num += 1
                    subtraj_first = subtraj_final
            else:
                return False
                    

    def _find_subtraj_final(self, traj, subtraj_first, ens_num):
        """
        Find the longest subtrajectory of trajectory which starts at
        subtraj_first and satifies self.ensembles[ens_num].can_append

        Returns
        -------
        int
            Frame of traj which is the final frame for a subtraj starting at
            subtraj_first and satisfying self.ensembles[ens_num]
        """
        subtraj_final = subtraj_first
        traj_final = len(traj)
        ens = self.ensembles[ens_num]
        subtraj = traj[slice(subtraj_first, subtraj_final+1)]
        # if we're in the ensemble or could eventually be in the ensemble,
        # we keep building the subtrajectory
        while ( (ens.can_append(subtraj) or ens(subtraj)) and 
                    subtraj_final < traj_final):
            subtraj_final += 1
            # TODO: replace with append; probably faster
            subtraj = traj[slice(subtraj_first, subtraj_final+1)]
        return subtraj_final
    
    def _find_subtraj_first(self, traj, subtraj_final, ens_num):
        subtraj_first = subtraj_final-1
        traj_first = 0
        ens = self.ensembles[ens_num]
        subtraj = traj[slice(subtraj_first, subtraj_final)]
        while ( (ens.can_prepend(subtraj) or ens(subtraj)) and
               subtraj_first >= traj_first):
            subtraj_first -= 1
            subtraj = traj[slice(subtraj_first, subtraj_final)]
        return subtraj_first+1


    def forward(self, trajectory):
        # treat this like we're implementing a regular expression parser ...
        # .*ensemble.+ ; but we have to do this for all possible matches
        # There are three tests we consider:
        # 1. subtraj_final - subtraj_first > 0: Do we obtain a subtrajectory?
        # 2. subtraj_final == traj_final: Have we assigned all the frames?
        # 3. ens_num == final_ens: are we looking at the last ensemble
        # Vaious combinations of these result in three possible outcomes:
        # (a) return True (we can append)
        # (b) return False (we can't append)
        # (c) loop around to text another subtrajectory (we can't tell)
        # Returning false can only happen if all ensembles have been tested
        traj_final = len(trajectory)
        final_ens = len(self.ensembles)-1
        #print traj_final, final_ens
        subtraj_first = 0
        ens_num = 0
        ens_first = 0
        while True: #  main loop, with various 
            subtraj_final = self._find_subtraj_final(trajectory, 
                                                     subtraj_first, ens_num)
            #print (ens_num,
                    #"("+str(subtraj_first)+","+str(subtraj_final)+")"
                  #)
            if subtraj_final - subtraj_first > 0:
                subtraj = trajectory[slice(subtraj_first, subtraj_final)]
                if ens_num == final_ens:
                    if subtraj_final == traj_final:
                        # we're in the last ensemble and the whole
                        # trajectory is assigned: can we append?
                        #print "Returning can_append"
                        return self.ensembles[ens_num].can_append(subtraj)
                    else:
                        #print "Returning false due to incomplete assigns:",
                        #print subtraj_final, "!=", traj_final
                        return False # in final ensemble, not all assigned
                else:
                    # subtraj existed, but not yet final ensemble
                    # so we start with the next ensemble
                    ens_num += 1
                    subtraj_first = subtraj_final
                    #print "Moving to the next ensemble", ens_num
            else:
                if subtraj_final == traj_final:
                    # all frames assigned, but not all ensembles finished;
                    # next frame might satisfy next ensemble
                    return True
                else:
                    # not all frames assigned, couldn't find a sequence
                    # start over with sequences that begin with the next
                    # ensemble
                    if ens_first == final_ens:
                        #print "Started with the last ensemble, got nothin'"
                        return False
                    else:
                        ens_first += 1
                        ens_num = ens_first
                        subtraj_first = 0
                    


    def backward(self, trajectory):
        # based on .forward(); see notes there for algorithm details
        traj_first = 0
        first_ens = 0
        subtraj_final = len(trajectory)
        ens_final = len(self.ensembles)-1
        ens_num = ens_final
        while True:
            subtraj_first = self._find_subtraj_first(trajectory,
                                                     subtraj_final, ens_num)
            #print (ens_num,
                    #"("+str(subtraj_first)+","+str(subtraj_final)+")"
                  #)
            if subtraj_final - subtraj_first > 0:
                subtraj = trajectory[slice(subtraj_first, subtraj_final)]
                if ens_num == first_ens:
                    if subtraj_first == traj_first:
                        return self.ensembles[ens_num].can_prepend(subtraj)
                    else:
                        return False
                else:
                    ens_num -= 1
                    subtraj_final = subtraj_first
            else:
                if subtraj_first == traj_first:
                    return True
                else:
                    if ens_final == first_ens:
                        return False
                    else:
                        ens_final -= 1
                        ens_num = ens_final
                        subtraj_final = len(trajectory)

    def __str__(self):
        head = "[\n"
        tail = "\n]"
        sequence_str = ",\n".join([str(ens) for ens in self.ensembles])
        return head+sequence_str+tail

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
    
    def __call__(self, trajectory, lazy=None):
        length = trajectory.frames
        if type(self.length) is int:
            return length == self.length
        else:
            return length >= self.length.start and (self.length.stop is None or length < self.length.stop)
        
    def forward(self, trajectory):
        length = trajectory.frames
        if type(self.length) is int:
            return length < self.length
        else:
            return self.length.stop is None or length < self.length.stop - 1

    def backward(self, trajectory):
        return self.forward(trajectory)

    can_append = forward
    can_prepend = backward

    def __str__(self):
        if type(self.length) is int:
            return 'len(x) = {0}'.format(self.length)
        else:
            return 'len(x) in [{0}, {1}]'.format(self.length.start, self.length.stop - 1)
        
class VolumeEnsemble(Ensemble):
    '''
    Describes an path ensemble using a volume object
    '''    
    def __init__(self, volume, frames = slice(1,-1), lazy = True):
        super(VolumeEnsemble, self).__init__()
        self._volume = volume
        self.frames = frames
        self.lazy = lazy
    
    @property
    def volume(self):
        '''
        The volume that is used in the specification.
        '''
        return self._volume
    
class InXEnsemble(VolumeEnsemble):
    '''
    Represents an ensemble where all the selected frames of the trajectory
    are in a specified volume
    '''

    def can_append(self, trajectory):
        return self(trajectory[slice(trajectory.frames-1, None)])

    def can_prepend(self, trajectory):
        return self(trajectory[slice(0,1)])
        
    
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
    
    def __call__(self, trajectory, lazy=None):
        if type(self.frames) is int:
            if trajectory.frames > self.frames and trajectory.frames >= -self.frames:
                return self.volume(trajectory[self.frames])            
            else:
                return True
        else:
            if (self.lazy and lazy is None) or (lazy and lazy is not None):
                for s in trajectory[self.frames][0::max(1,trajectory.frames - 1)]:
                    if not self.volume(s):
                        return False
            else:
                for s in trajectory[self.frames]:
                    if not self.volume(s):
                        return False
                        
            return True

    def __str__(self):
        if type(self.frames) is int:
            return 'x[{0}] in {1}'.format(str(self.frames), str(self.volume))
        else:
            return 'x[t] in {2} for all t in [{0}:{1}])'.format(self.frames.start, self.frames.stop, self.volume)


class OutXEnsemble(InXEnsemble):
    '''
    Represents an ensemble where all the selected frames from the trajectory
    are outside a specified volume
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
    Represents an ensemble where at least one of the selected frames from
    the trajectory visit a specified volume
    '''

    def __str__(self):
        if type(self.frames) is int:
            return 'x[{0}] in {1}'.format(str(self.frames), str(self.volume))
        else:
            return 'x[t] in {2} for one t in [{0}:{1}]'.format(self.frames.start, self.frames.stop, self.volume)

    def can_append(self, trajectory):
        return True

    def can_prepend(self, trajectory):
        return True


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

    
    def __call__(self, trajectory, lazy=None):
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
    Represents an ensemble where at least one frame of the trajectory is
    outside the specified volume
    '''
    def __str__(self):
        if type(self.frames) is int:
            return 'x[{0}] not in {1}'.format(str(self.frames), str(self.volume))
        else:
            return 'x[t] in {2} for one t in [{0}:{1}])'.format(self.frames.start, self.frames.stop, self.volume)
        
    @property
    def volume(self):
        return ~ self._volume


class ExitsXEnsemble(VolumeEnsemble):
    """
    Represents an ensemble where two successive frames from the selected
    frames of the trajectory crossing from inside to outside the given volume.
    """
    def __init__(self, volume, frames = slice(None), lazy=False):
        # changing the defaults for frames and lazy; prevent single frame
        if type(frames) is int:
            raise ValueError('Exits/EntersXEnsemble require more than one frame')
        super(ExitsXEnsemble, self).__init__(volume,frames,lazy)

    def __str__(self):
        domain = 'exists x[t], x[t+1] in [{0}:{1}] '.format(
                            self.frames.start, self.frames.stop )
        result = 'such that x[t] in {0} and x[t+1] not in {0}'.format(
                            self.volume)
        return domain+result

    def __call__(self, trajectory, lazy=None):
        if type(self.frames) is int:
            # in case you changed self.frames after intialization
            raise ValueError('ExitsXEnsemble requires more than one frame')
        else:
            subtraj = trajectory[self.frames]
            for i in range(len(subtraj)-1):
                frame_i = subtraj[i]
                frame_iplus = subtraj[i+1]
                if self.volume(frame_i) and not self.volume(frame_iplus):
                    return True
        return False

class EntersXEnsemble(ExitsXEnsemble):
    """
    Represents an ensemble where two successive frames from the selected
    frames of the trajectory crossing from outside to inside the given volume.
    """
    def __str__(self):
        domain = 'exists x[t], x[t+1] in [{0}:{1}] '.format(
                            self.frames.start, self.frames.stop )
        result = 'such that x[t] not in {0} and x[t+1] in {0}'.format(
                            self.volume)
        return domain+result

    def __call__(self, trajectory, lazy=None):
        if type(self.frames) is int:
            # in case you changed self.frames after intialization
            raise ValueError('EntersXEnsemble requires more than one frame')
        else:
            subtraj = trajectory[self.frames]
            for i in range(len(subtraj)-1):
                frame_i = subtraj[i]
                frame_iplus = subtraj[i+1]
                if not self.volume(frame_i) and self.volume(frame_iplus):
                    return True
        return False
            
        
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
        
    def __call__(self, trajectory, lazy=None):
        return self.ensemble(self._alter(trajectory), lazy)

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
        return (LengthEnsemble(slice(3,None)) & InXEnsemble(volume_a, 0) & InXEnsemble(volume_b, -1)) & OutXEnsemble(volume_a | volume_b, slice(1,-1), lazy)

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

        return (LengthEnsemble(slice(3,None)) & InXEnsemble(volume_a, 0) & InXEnsemble(volume_b, -1)) & (LeaveXEnsemble(volume_x) & OutXEnsemble(volume_a | volume_b, slice(1,-1), lazy))
    
    hbar = StartXEnsemble
