'''
Created on 03.09.2014

@author: jan-hendrikprinz, David W.H. Swenson
'''

from openpathsampling.todict import restores_as_full_object

import openpathsampling as paths

import logging
from ops_logging import initialization_logging
logger = logging.getLogger(__name__)
init_log = logging.getLogger('opentis.initialization')

# TODO: Make Full and Empty be Singletons to avoid storing them several times!

@restores_as_full_object
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
    Maybe replace - by / to get better notation. So far it has not been used
    '''

    use_shortcircuit = True

    def __init__(self):
        '''
        A path volume defines a set of paths.
        '''

#        self._traj = dict()
#        self.last = None
        self.name = ''

    def __eq__(self, other):
        if self is other:
            return True
        return str(self) == str(other)

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

    def oom_matrix(self, oom):
        """
        Return the oom representation where the OOM is based on a set of volumes

        """

        # Needs to be implemented by the actual class

        return None
    
    def can_append(self, trajectory):
        '''
        Returns true, if the trajectory so far can still be in the ensemble
        if it is appended by a frame. To check, it assumes that the
        trajectory to length L-1 is okay. This is mainly for interactive
        usage, when a trajectory is generated.
        
        Parameters
        ----------
        trajectory : Trajectory
            the actual trajectory to be tested
        
        Returns
        -------
        can_append : bool
            Returns true or false if using a forward step (extending the
            trajectory forward in time at its end) `trajectory` could  still
            be in the ensemble and thus makes sense to continue a simulation
        

        Notes
        -----
        This is only tricky for this that depend on the history like
        HitXEnsemble or LeaveXEnsembles. In theory these can only be checked
        if the full range of frames has been generated. This could be
        triggered, when the last frame is reached.  This is even more
        difficult if this depends on the length.
        '''
        return True        
    
    def can_prepend(self, trajectory):
        '''
        Returns true, if the trajectory so far can still be in the ensemble
        if it is prepended by a frame. To check, it assumes that the
        trajectory from index 1 is okay. This is mainly for interactive
        usage, when a trajectory is generated using a backward move.
        
        Parameters
        ----------
        trajectory : Trajectory
            the actual trajectory to be tested
        
        Returns
        -------
        can_prepend : bool
            Returns true or false if using a backward step (extending the
            trajectory backwards in time at its beginning) `trajectory`
            could  still be in the ensemble and thus makes sense to continue
            a simulation
        
        Notes
        
        This is only tricky for this that depend on the history like
        HitXEnsemble or LeaveXEnsembles. In theory these can only be checked
        if the full range of frames has been generated. This could be
        triggered, when the last frame is reached.  This is even more
        difficult if this depends on the length.
        '''
        return True        



    def find_valid_slices(self, trajectory, lazy=True, 
                          max_length=None, min_length=1, overlap=1):
        '''
        Returns a list of trajectories that contain sub-trajectories which
        are in the given ensemble.

        Parameters
        ----------
        trajectory : Trajectory
            the actual trajectory to be splitted into ensemble parts
        lazy : boolean
            if True will use a faster almost linear algorithm, while False
            will run through all possibilities starting with the largest
            ones
        max_length : int > 0
            if set this determines the maximal size to be tested (is mainly
            used in the recursion)
        min_length : int > 0
            if set this determines the minimal size to be tested (in lazy
            mode might no
        overlap : int >= 0
            determines the allowed overlap of all trajectories to be found.
            A value of x means that two sub-trajectorie can share up to x
            frames at the beginning and x frames at the end.  Default is 1

        Returns
        -------
        list of slices
            Returns a list of index-slices for sub-trajectories in
            trajectory that are in the ensemble.
        '''
        ensemble_list = []

        length = len(trajectory)

        if max_length is None:
            max_length = length

        max_length = min(length, max_length)
        min_length = max(1, min_length)

        logger.debug("Looking for subtrajectories in " + str(trajectory))

        if not lazy:
            # this tries all possible sub-trajectories starting with the
            # longest ones and uses recursion
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
                            list_left = self.find_valid_slices(tt_left, 
                                                               max_length=l)

                            tt_right = trajectory[start + l - pad:length]
                            list_right = self.find_valid_slices(tt_right, 
                                                                max_length=l)

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
                if self.can_append(tt) and end<length:
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
        Returns a list of trajectories that contain sub-trajectories which
        are in the given ensemble.

        Parameters
        ----------
        trajectory : Trajectory
            the actual trajectory to be splitted into ensemble parts
        lazy : boolean
            if True will use a faster almost linear algorithm, while False
            will run through all possibilities starting with the largest
            ones
        max_length : int > 0
            if set this determines the maximal size to be tested (is mainly
            used in the recursion)
        min_length : int > 0
            if set this determines the minimal size to be tested (in lazy
            mode might no
        overlap : int >= 0
            determines the allowed overlap of all trajectories to be found.
            A value of x means that two sub-trajectory can share up to x
            frames at the beginning and x frames at the end.  Default is 1

        Returns
        -------
        list of Trajectory
            Returns a list of sub-trajectories in trajectory that are in the
            ensemble.

        Notes
        -----
        This uses self.find_valid_slices and returns the actual sub-trajectories
        '''

        indices = self.find_valid_slices(trajectory, lazy, max_length, 
                                         min_length, overlap)

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
        elif type(other) is EmptyEnsemble:
            return self
        elif type(other) is FullEnsemble:
            return other
        else:
            return OrEnsemble(self, other)

    def __xor__(self, other):
        if self is other:
            return EmptyEnsemble()
        elif type(other) is EmptyEnsemble:
            return self
        elif type(other) is FullEnsemble:
            return NegatedEnsemble(self)        
        else:
            return XorEnsemble(self, other)

    def __and__(self, other):
        if self is other:
            return self
        elif type(other) is EmptyEnsemble:
            return other
        elif type(other) is FullEnsemble:
            return self
        else:
            return AndEnsemble(self, other)

    def __sub__(self, other):
        if self is other:
            return EmptyEnsemble()
        elif type(other) is EmptyEnsemble:
            return self
        elif type(other) is FullEnsemble:
            return EmptyEnsemble()
        else:
            return SubEnsemble(self, other)
        
    def __invert__(self):
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
                
class LoadedEnsemble(Ensemble):
    '''
    Represents an ensemble the contains trajectories of a specific length
    '''
    def __init__(self, name, description):
        '''
        A path ensemble that describes path of a specific length

        Parameters
        ----------
        length : int or slice
            The specific length (int) or the range of allowed trajectory lengths (slice)
        '''

        super(LoadedEnsemble, self).__init__()

        self.name = name
        self.description = description
        pass

    def __str__(self):
        return self.description

@restores_as_full_object
class EmptyEnsemble(Ensemble):
    '''
    The empty path ensemble of no trajectories.
    '''
    def __init__(self):
        super(EmptyEnsemble, self).__init__()

    def __call__(self, trajectory, lazy=None):
        return False

    def can_append(self, trajectory):
        return False

    def can_prepend(self, trajectory):
        return False

    def __invert__(self):
        return FullEnsemble()

    def __sub__(self, other):
        return EmptyEnsemble()

    def __and__(self, other):
        return self

    def __xor__(self, other):
        return other

    def __or__(self, other):
        return other

    def __str__(self):
        return 'empty'

    def oom_matrix(self, oom):
        # Zero matrix
        return None

@restores_as_full_object
class FullEnsemble(Ensemble):
    '''
    The full path ensemble of all possible trajectories.
    '''
    def __init__(self):
        super(FullEnsemble, self).__init__()

    def __call__(self, trajectory, lazy=None):
        return True
    
    def can_append(self, trajectory):
        return True

    def can_prepend(self, trajectory):
        return True

    def __invert__(self):
        return EmptyEnsemble()

    def __sub__(self, other):
        if type(other) is EmptyEnsemble:
            return self
        elif type(other) is FullEnsemble:
            return EmptyEnsemble()
        else:
            return NegatedEnsemble(other)

    def __and__(self, other):
        return other

    def __xor__(self, other):
        if type(other) is EmptyEnsemble:
            return self
        elif type(other) is FullEnsemble:
            return EmptyEnsemble()
        else:
            return NegatedEnsemble(other)

    def __or__(self, other):
        return self

    def __str__(self):
        return 'all'

    def oom_matrix(self, oom):
        # Full matrix
        return None

@restores_as_full_object
class NegatedEnsemble(Ensemble):
    '''
    Negates an Ensemble and simulates a `not` statement
    '''
    def __init__(self, volume):
        super(NegatedEnsemble, self).__init__()
        self.ensemble = volume
        
    def __call__(self, trajectory, lazy=None):
        return not self.ensemble(trajectory, lazy)

    def can_append(self, trajectory):
        # We cannot guess the result here so keep on running forever
        return True

    def can_prepend(self, trajectory):
        # We cannot guess the result here so keep on running forever
        return True

    def __str__(self):
        return 'not ' + str(self.ensemble2)

@restores_as_full_object
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

    def to_dict(self):
        return { 'ensemble1' : self.ensemble1, 'ensemble2' : self.ensemble2 }

    def __call__(self, trajectory, lazy=None):
        # Shortcircuit will automatically skip the second part of the combination if the result does not depend on it!
        # This makes sense since the expensive part is the ensemble testing not computing two logic operations
        if Ensemble.use_shortcircuit:
            a = self.ensemble1(trajectory, lazy)
            logger.debug("Combination: " + self.ensemble1.__class__.__name__ + 
                         " is "+str(a))
            logger.debug("Combination: " + self.ensemble2.__class__.__name__ +   
                         " is " +str(self.ensemble2(trajectory, lazy)))
            logger.debug("Combination: returning " + 
                         str(self.fnc(a,self.ensemble2(trajectory,lazy))))
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

    # Forward / Backward is tricky
    # We can do the following. If a or b is true this means that the real result could be false or true, we just
    # keep going but we should have stopped. If a or b is false this means for that ensemble continuing is not
    # feasible and so false really means false. To check if a logical combination should be continued just
    # try for all true values a potential false and check if we should continue.

    def _continue_fnc(self, a, b):
        fnc = self.fnc
        res = fnc(a,b)
        if a is True:
            res |= fnc(False, b)
        if b is True:
            res |= fnc(a, False)
        if a is True and b is True:
            res |= fnc(False, False)

        return res

    def can_append(self, trajectory):
        if Ensemble.use_shortcircuit:
            a = self.ensemble1.can_append(trajectory)
            res_true = self._continue_fnc(a, True)
            res_false = self._continue_fnc(a, False)
            if res_false == res_true:
                # result is independent of ensemble_b so ignore it
                return res_true
            else:
                b = self.ensemble2.can_append(trajectory)
                if b is True:
                    return res_true
                else:
                    return res_false
        else:
            return self.fnc(self.ensemble1.can_append(trajectory), self.ensemble2.can_append(trajectory))

    def can_prepend(self, trajectory):
        if Ensemble.use_shortcircuit:
            a = self.ensemble1.can_prepend(trajectory)
            res_true = self._continue_fnc(a, True)
            res_false = self._continue_fnc(a, False)
            if res_false == res_true:
                # result is independent of ensemble_b so ignore it
                return res_true
            else:
                b = self.ensemble2.can_prepend(trajectory)
                if b is True:
                    return res_true
                else:
                    return res_false
        else:
            return self.fnc(self.ensemble1.can_prepend(trajectory), self.ensemble2.can_prepend(trajectory))

    def __str__(self):
#        print self.sfnc, self.ensemble1, self.ensemble2, self.sfnc.format('(' + str(self.ensemble1) + ')' , '(' + str(self.ensemble1) + ')')
        return self.sfnc.format('(\n' + Ensemble._indent(str(self.ensemble1)) + '\n)' , '(\n' + Ensemble._indent(str(self.ensemble2)) + '\n)')

@restores_as_full_object
class OrEnsemble(EnsembleCombination):
    def __init__(self, ensemble1, ensemble2):
        super(OrEnsemble, self).__init__(ensemble1, ensemble2, fnc = lambda a,b : a or b, str_fnc = '{0}\nor\n{1}')

@restores_as_full_object
class AndEnsemble(EnsembleCombination):
    def __init__(self, ensemble1, ensemble2):
        super(AndEnsemble, self).__init__(ensemble1, ensemble2, fnc = lambda a,b : a and b, str_fnc = '{0}\nand\n{1}')

@restores_as_full_object
class XorEnsemble(EnsembleCombination):
    def __init__(self, ensemble1, ensemble2):
        super(XorEnsemble, self).__init__(ensemble1, ensemble2, fnc = lambda a,b : a ^ b, str_fnc = '{0}\nxor\n{1}')

@restores_as_full_object
class SubEnsemble(EnsembleCombination):
    def __init__(self, ensemble1, ensemble2):
        super(SubEnsemble, self).__init__(ensemble1, ensemble2, fnc = lambda a,b : a and not b, str_fnc = '{0}\nand not\n{1}')

@restores_as_full_object
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

    Notes
    -----
        TODO: Overlap features not implemented because ohmygod this was hard
        enough already.
    """

    def __init__(self, ensembles, min_overlap=0, max_overlap=0, greedy=False):
        # make tuples of the min/max overlaps
        super(SequentialEnsemble, self).__init__()
        if type(min_overlap) is int:
            min_overlap = (min_overlap, )*(len(ensembles)-1)
        if type(max_overlap) is int:
            max_overlap = (max_overlap, )*(len(ensembles)-1)

        self.ensembles = ensembles
        self.min_overlap = min_overlap
        self.max_overlap = max_overlap
        self.greedy = greedy

        # sanity checks
        if len(self.min_overlap) != len(self.max_overlap):
            raise ValueError("len(min_overlap) != len(max_overlap)")
        if len(self.min_overlap) != len(self.ensembles)-1:
            raise ValueError("Number of overlaps doesn't match number of transitions")
        for i in range(len(self.min_overlap)):
            if min_overlap[i] > max_overlap[i]:
                raise ValueError("min_overlap greater than max_overlap!")

    def transition_frames(self, trajectory, lazy=None):
        # it is easiest to understand this decision tree as a simplified
        # version of the can_append decision tree; see that for detailed
        # comments
        ens_num = 0
        subtraj_first = 0
        traj_final = len(trajectory)
        final_ens = len(self.ensembles)-1
        transitions = []
        while True:
            if ens_num <= final_ens:
                subtraj_final = self._find_subtraj_final(trajectory, 
                                                         subtraj_first, ens_num)
            else:
                return transitions
            if subtraj_final - subtraj_first > 0:
                subtraj = trajectory[slice(subtraj_first, subtraj_final)]
                if ens_num == final_ens:
                    if subtraj_final == traj_final:
                        # success
                        transitions.append(subtraj_final)
                        return transitions
                    else:
                        # fails because we have more frames to assign
                        transitions.append(subtraj_final) 
                        return transitions
                else:
                    ens_num += 1
                    transitions.append(subtraj_final)
                    subtraj_first = subtraj_final
            else:
                if ens_num <= final_ens and self.ensembles[ens_num](paths.Trajectory([])):
                    ens_num += 1
                    transitions.append(subtraj_final)
                    subtraj_first = subtraj_final
                else:
                    return transitions


    def __call__(self, trajectory, lazy=None):
        logger.debug("Looking for transitions in trajectory " + str(trajectory))
        transitions = self.transition_frames(trajectory, lazy)
        logger.debug("Found transitions: " + str(transitions))
        # if we don't have the right number of transitions, or if the last 
        #print transitions
        if len(transitions) != len(self.ensembles):
            #print "Returns false b/c not enough ensembles"
            return False
        elif transitions[-1] != len(trajectory):
            #print "Returns false b/c not all frames assigned"
            return False

        subtraj_first = 0
        subtraj_i = 0
        while subtraj_i < len(self.ensembles):
            subtraj_final = transitions[subtraj_i]
            subtraj = trajectory[slice(subtraj_first, subtraj_final)]
            if self.ensembles[subtraj_i](subtraj) == False:
                #print "Returns false b/c ensemble",subtraj_i," fails"
                return False
            subtraj_i += 1
            subtraj_first = subtraj_final
        return True

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


    def can_append(self, trajectory):
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

        # logging startup
        logger.debug("Beginning can_append")
        for ens in self.ensembles:
            logger.debug(
                "Ensemble " + str(self.ensembles.index(ens)) + 
                " : " + ens.__class__.__name__
            )

        while True: #  main loop, with various 
            subtraj_final = self._find_subtraj_final(trajectory, 
                                                     subtraj_first, ens_num)
            logger.debug(
                str(ens_num) + " : " +
                "("+str(subtraj_first)+","+str(subtraj_final)+")"
            )
            if subtraj_final - subtraj_first > 0:
                subtraj = trajectory[slice(subtraj_first, subtraj_final)]
                if ens_num == final_ens:
                    if subtraj_final == traj_final:
                        # we're in the last ensemble and the whole
                        # trajectory is assigned: can we append?
                        logger.debug("Returning can_append for " + str(self.ensembles[ens_num].__class__.__name__))
                        return self.ensembles[ens_num].can_append(subtraj)
                    else:
                        logger.debug(
                            "Returning false due to incomplete assigns: " + 
                            str(subtraj_final) + "!=" + str(traj_final)
                        )
                        return False # in final ensemble, not all assigned
                else:
                    # subtraj existed, but not yet final ensemble
                    # so we start with the next ensemble
                    if subtraj_final != traj_final and not self.ensembles[ens_num](subtraj):
                        logger.debug(
                            "Couldn't assign frames " + str(subtraj_first) +
                            " through " + str(subtraj_final) + 
                            " to ensemble " + str(ens_num) + ": No match"
                        )

                    logger.debug(
                        "Assigning frames " + str(subtraj_first) +
                        " through " + str(subtraj_final) + 
                        " to ensemble " + str(ens_num)
                    )
                    ens_num += 1
                    subtraj_first = subtraj_final

                    logger.debug("Moving to the next ensemble " + str(ens_num))
            else:
                if subtraj_final == traj_final:
                    # all frames assigned, but not all ensembles finished;
                    # next frame might satisfy next ensemble
                    logger.debug("All frames assigned, more ensembles to go: returning True")
                    return True
                elif self.ensembles[ens_num](paths.Trajectory([])):
                    logger.debug("Moving on because of allowed zero-length ensemble")
                    ens_num += 1
                    subtraj_first = subtraj_final
                else:
                    # not all frames assigned, couldn't find a sequence
                    # start over with sequences that begin with the next
                    # ensemble
                    if ens_first == final_ens:
                        logger.debug("Started with the last ensemble, got nothin'")
                        return False
                    else:
                        logger.debug(
                            "Reassigning all frames, starting with ensemble " +
                            str(ens_first)
                        )
                        ens_first += 1
                        ens_num = ens_first
                        subtraj_first = 0


    def can_prepend(self, trajectory):
        # based on .can_append(); see notes there for algorithm details
        traj_first = 0
        first_ens = 0
        subtraj_final = len(trajectory)
        ens_final = len(self.ensembles)-1
        ens_num = ens_final

        # logging startup
        logger.debug("Beginning can_prepend")
        for i in range(len(self.ensembles)):
            logger.debug(
                "Ensemble " + str(i) + 
                " : " + self.ensembles[i].__class__.__name__
            )

        while True:
            subtraj_first = self._find_subtraj_first(trajectory,
                                                     subtraj_final, ens_num)
            logger.debug(
                str(ens_num) + " : " +
                "("+str(subtraj_first)+","+str(subtraj_final)+")"
            )
            if subtraj_final - subtraj_first > 0:
                subtraj = trajectory[slice(subtraj_first, subtraj_final)]
                if ens_num == first_ens:
                    if subtraj_first == traj_first:
                        logger.debug("Returning can_prepend")
                        return self.ensembles[ens_num].can_prepend(subtraj)
                    else:
                        logger.debug(
                            "Returning false due to incomplete assigns: " + 
                            str(subtraj_first) + "!=" + str(traj_first)
                        )
                        return False
                else:
                    if subtraj_first != traj_first and not self.ensembles[ens_num](subtraj):
                        logger.debug(
                            "Couldn't assign frames " + str(subtraj_first) +
                            " through " + str(subtraj_final) + 
                            " to ensemble " + str(ens_num) + ": No match"
                        )
                    logger.debug(
                        "Assigning frames " + str(subtraj_first) +
                        " through " + str(subtraj_final) + 
                        " to ensemble " + str(ens_num)
                    )
                    ens_num -= 1
                    subtraj_final = subtraj_first
                    logger.debug("Moving to the next ensemble " + str(ens_num))
            else:
                if subtraj_first == traj_first:
                    logger.debug("All frames assigned, more ensembles to go: returning True")
                    return True
                elif self.ensembles[ens_num](paths.Trajectory([])):
                    logger.debug("Moving on because of allowed zero-length ensemble")
                    ens_num -= 1
                    subtraj_final = subtraj_first
                else:
                    if ens_final == first_ens:
                        logger.debug("Started with the last ensemble, got nothin'")
                        return False
                    else:
                        logger.debug(
                            "Reassigning all frames, starting with ensemble " +
                            str(ens_final)
                        )
                        ens_final -= 1
                        ens_num = ens_final
                        subtraj_final = len(trajectory)

    def __str__(self):
        head = "[\n"
        tail = "\n]"
        sequence_str = ",\n".join([str(ens) for ens in self.ensembles])
        return head+sequence_str+tail


@restores_as_full_object
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

        super(LengthEnsemble, self).__init__()
        self.length = length
        pass
    
    def __call__(self, trajectory, lazy=None):
        length = trajectory.frames
        if type(self.length) is int:
            return length == self.length
        else:
            return length >= self.length.start and (self.length.stop is None or length < self.length.stop)
        
    def can_append(self, trajectory):
        length = trajectory.frames
        if type(self.length) is int:
            return length < self.length
        else:
            return self.length.stop is None or length < self.length.stop - 1

    def can_prepend(self, trajectory):
        return self.can_append(trajectory)

    def __str__(self):
        if type(self.length) is int:
            return 'len(x) = {0}'.format(self.length)
        else:
            start = self.length.start
            if start is None:
                start = 0
            stop = self.length.stop
            if stop is None:
                stop = 'infty'
            else:
                stop = str(self.length.stop - 1)
            return 'len(x) in [{0}, {1}]'.format(start, stop)

@restores_as_full_object
class VolumeEnsemble(Ensemble):
    '''
    Describes an path ensemble using a volume object
    '''    
    def __init__(self, volume, lazy = True):
        super(VolumeEnsemble, self).__init__()
        self.volume = volume
        self.lazy = lazy

    @property
    def _volume(self):
        '''
        The volume that is used in the specification.
        '''
        return self.volume

@restores_as_full_object
class InXEnsemble(VolumeEnsemble):
    '''
    Represents an ensemble where all the selected frames of the trajectory
    are in a specified volume
    '''

    def can_append(self, trajectory):
        if len(trajectory) == 0:
            return True
        else:
            return self(trajectory[slice(trajectory.frames-1, None)])

    def can_prepend(self, trajectory):
        if len(trajectory) == 0:
            return True
        else:
            return self(trajectory[slice(0,1)])
        
    
    def __call__(self, trajectory, lazy=None):
        if len(trajectory) == 0:
            return False
        for frame in trajectory:
            if not self._volume(frame):
                return False
        return True

    def __invert__(self):
        return LeaveXEnsemble(self.volume, self.frames, self.lazy)

    def __str__(self):
        return 'x[t] in {0} for all t'.format(self._volume)


@restores_as_full_object
class OutXEnsemble(InXEnsemble):
    '''
    Represents an ensemble where all the selected frames from the trajectory
    are outside a specified volume
    '''    
    @property
    def _volume(self):
        return ~ self.volume
    
    def __str__(self):
        return 'x[t] in {0} for all t'.format(self._volume)

    def __invert__(self):
        return HitXEnsemble(self.volume, self.frames, self.lazy)

@restores_as_full_object
class HitXEnsemble(VolumeEnsemble):
    '''
    Represents an ensemble where at least one of the selected frames from
    the trajectory visit a specified volume
    '''

    def __str__(self):
        return 'exists t such that x[t] in {0}'.format(self._volume)

    def __call__(self, trajectory, lazy=None):
        '''
        Returns True if the trajectory is part of the PathEnsemble
        
        Parameters
        ----------
        trajectory : Trajectory
            The trajectory to be checked
        '''
        for frame in trajectory:
            if self._volume(frame):
                return True
        return False

    def __invert__(self):
        return OutXEnsemble(self.volume, self.frames, self.lazy)

@restores_as_full_object
class LeaveXEnsemble(HitXEnsemble):
    '''
    Represents an ensemble where at least one frame of the trajectory is
    outside the specified volume
    '''
    def __str__(self):
        return 'exists t such that x[t] in {0}'.format(self._volume)
      
    @property
    def _volume(self):
        # effectively use HitXEnsemble but with inverted volume
        return ~ self.volume

    def __invert__(self):
        return InXEnsemble(self.volume, self.frames, self.lazy)

    def __call__(self, trajectory, lazy=None):
        for frame in trajectory:
            if self._volume(frame):
                return True
        return False

@restores_as_full_object
class ExitsXEnsemble(VolumeEnsemble):
    """
    Represents an ensemble where two successive frames from the selected
    frames of the trajectory crossing from inside to outside the given volume.
    """
    def __init__(self, volume, lazy=False):
        # changing the defaults for frames and lazy; prevent single frame
        super(ExitsXEnsemble, self).__init__(volume, lazy)

    def __str__(self):
        domain = 'exists x[t], x[t+1] '
        result = 'such that x[t] in {0} and x[t+1] not in {0}'.format(
                            self._volume)
        return domain+result

    def __call__(self, trajectory, lazy=None):
        subtraj = trajectory
        for i in range(len(subtraj)-1):
            frame_i = subtraj[i]
            frame_iplus = subtraj[i+1]
            if self._volume(frame_i) and not self._volume(frame_iplus):
                return True
        return False

@restores_as_full_object
class EntersXEnsemble(ExitsXEnsemble):
    """
    Represents an ensemble where two successive frames from the selected
    frames of the trajectory crossing from outside to inside the given volume.
    """
    def __str__(self):
        domain = 'exists x[t], x[t+1] '
        result = 'such that x[t] not in {0} and x[t+1] in {0}'.format(
                            self._volume)
        return domain+result

    def __call__(self, trajectory, lazy=None):
        subtraj = trajectory
        for i in range(len(subtraj)-1):
            frame_i = subtraj[i]
            frame_iplus = subtraj[i+1]
            if not self._volume(frame_i) and self._volume(frame_iplus):
                return True
        return False

@restores_as_full_object
class WrappedEnsemble(Ensemble):
    '''
    Represents an ensemble where an altered version of a trajectory (extended, reversed, cropped) is part of a given ensemble
    '''
    def __init__(self, ensemble):
        '''
        Represents an ensemble which is the given ensemble but for trajectories where some trajectory is prepended
        '''
        super(WrappedEnsemble, self).__init__()
        self.ensemble = ensemble

        # you can also build wrapped ensembles with more flexibility when using
        # a property for _new_ensemble
        self._new_ensemble = self.ensemble

    def __call__(self, trajectory, lazy=None):
        return self._new_ensemble(self._alter(trajectory), lazy)

    def _alter(self, trajectory):
        return trajectory
        
    def can_append(self, trajectory):
        return self._new_ensemble.can_append(self._alter(trajectory))

    def can_prepend(self, trajectory):
        return self._new_ensemble.can_prepend(self._alter(trajectory))

@restores_as_full_object
class SlicedTrajectoryEnsemble(WrappedEnsemble):
    '''
    An ensemble which alters the trajectory by looking at a given Python
    slice of the list of frames.
    '''
    def __init__(self, ensemble, aslice):
        super(SlicedTrajectoryEnsemble, self).__init__(ensemble)
        if type(aslice) == int:
            if aslice == -1:
                self.slice = slice(aslice,None)
            else:
                self.slice = slice(aslice, aslice+1)
        else:
            self.slice = aslice

    def _alter(self, trajectory):
        return trajectory[self.slice]

    def __str__(self):
        # TODO: someday may add different string support for slices with
        # only one frame
        start = "" if self.slice.start is None else str(self.slice.start)
        stop = "" if self.slice.stop is None else str(self.slice.stop)
        step = "" if self.slice.step is None else " every "+str(self.slice.step)
        return ("(" + self.ensemble.__str__() +
                " in {" + start + ":" + stop + "}" + step + ")")


@restores_as_full_object
class BackwardPrependedTrajectoryEnsemble(WrappedEnsemble):
    '''
    Represents an ensemble which is the given ensemble but for trajectories where some trajectory is prepended
    '''
    def __init__(self, ensemble, trajectory):        
        super(BackwardPrependedTrajectoryEnsemble, self).__init__(ensemble)
        self.add_traj = trajectory        

    def _alter(self, trajectory):
#        print [ s.idx for s in trajectory.reversed + self.add_traj]
        return trajectory.reversed + self.add_traj

@restores_as_full_object
class ForwardAppendedTrajectoryEnsemble(WrappedEnsemble):
    '''
    Represents an ensemble which is the given ensemble but for trajectories where some trajectory is appended
    '''
    def __init__(self, ensemble, trajectory):
        super(ForwardAppendedTrajectoryEnsemble, self).__init__(ensemble)
        self.add_traj = trajectory

    def _alter(self, trajectory):
        return self.add_traj + trajectory

@restores_as_full_object
class ReversedTrajectoryEnsemble(WrappedEnsemble):
    '''
    Represents an ensemble 
    '''
    def _alter(self, trajectory):
        return trajectory.reverse()

@restores_as_full_object
class AppendedNameEnsemble(WrappedEnsemble):
    '''
    Adds string to ensemble name: necessary to have multiple copies of an ensemble.
    '''
    def __init__(self, ensemble, label):
        self._label = label
        super(AppendedNameEnsemble, self).__init(ensemble)

    def __str__(self):
        return self.ensemble.__str__() + " " + self.label



@restores_as_full_object
class OptionalEnsemble(WrappedEnsemble):
    '''
    Makes it optional to satisfy a given ensemble (primarily useful in
    SequentialEnsembles)
    '''

    def __init__(self, ensemble):
        super(OptionalEnsemble, self).__init__(ensemble)
        self._new_ensemble = LengthEnsemble(0) | self.ensemble

    def __str__(self):
        return "{"+self.ensemble.__str__()+"} (OPTIONAL)"

@restores_as_full_object
class SingleFrameEnsemble(WrappedEnsemble):
    '''
    Convenience ensemble to `and` a LengthEnsemble(1) with a given ensemble.
    Frequently used for SequentialEnsembles.

    Attributes
    ----------
    ensemble : Ensemble
        the ensemble which should be represented in the single frame

    Notes
    -----
    We allow the user to choose to be stupid: if, for example, the user
    tries to make a SingleFrameEnsemble from an ensemble which requires
    more than one frame to be satisfied (e.g., a SequentialEnsemble with
    more than one subensemble), it can be created, but no path will ever
    satisfy it. Since we can't stop all possible mistakes, we don't bother
    here.
    '''
    def __init__(self, ensemble):
        super(SingleFrameEnsemble, self).__init__(ensemble)
        self._new_ensemble = LengthEnsemble(1) & self.ensemble

    def __str__(self):
        return "{"+self.ensemble.__str__()+"} (SINGLE FRAME)"

@restores_as_full_object
class MinusInterfaceEnsemble(SequentialEnsemble):
    '''
    This creates an ensemble for the minus interface. 

    Parameters
    ----------
    state_vol : Volume
        The Volume which defines the state for this minus interface
    innermost_vol : Volume
        The Volume defining the innermost interface with which this minus
        interface does its replica exchange.
    n_l : integer (greater than one)
        The number of segments crossing innermost_vol for this interface.
    
    The specific implementation allows us to use the multiple-segment minus
    ensemble described by Swenson and Bolhuis. The minus interface was
    originally developed by van Erp. For more details, see the section
    "Anatomy of a PathMover: the Minus Move" in the OpenPathSampling
    Documentation.

    References
    ----------
    D.W.H. Swenson and P.G. Bolhuis. J. Chem. Phys. 141, 044101 (2014). 
    doi:10.1063/1.4890037
    '''

    # don't store unnecessary stuff we recreate at initialization
    # TODO: Check with David if it makes sense to store these and allow
    # them being used in __init__ instead of the self-made ones

    _excluded_attr = ['ensembles', 'min_overlap', 'max_overlap']

    def __init__(self, state_vol, innermost_vol, n_l=2, greedy=False):
        if (n_l < 2):
            raise ValueError("The number of segments n_l must be at least 2")

        self.state_vol = state_vol
        self.innermost_vol = innermost_vol
        self.greedy = greedy
        inA = InXEnsemble(state_vol)
        outA = OutXEnsemble(state_vol)
        outX = OutXEnsemble(innermost_vol)
        inX = InXEnsemble(innermost_vol)
        leaveX = LeaveXEnsemble(innermost_vol)
        interstitial = outA & inX
        self._segment_ensemble = EnsembleFactory.TISEnsemble(
            state_vol, state_vol, innermost_vol)
        #interstitial = InXEnsemble(innermost_vol - state_vol)
        start = [
            SingleFrameEnsemble(inA),
            OptionalEnsemble(interstitial),
        ]
        loop = [
            outA & leaveX,
            inX # & hitA # redundant due to stop req for previous outA
        ]
        end = [
            outA & leaveX,
            OptionalEnsemble(interstitial),
            SingleFrameEnsemble(inA)
        ]
        ensembles = start + loop*(n_l-1) + end

        self.n_l = n_l

        super(MinusInterfaceEnsemble, self).__init__(ensembles, greedy=greedy)

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
        # TODO: this is actually only for flexible path length TPS now
        return SequentialEnsemble([
            SingleFrameEnsemble(InXEnsemble(volume_a)),
            OutXEnsemble(volume_a | volume_b),
            SingleFrameEnsemble(InXEnsemble(volume_b))
        ])



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
        ens = SequentialEnsemble([
            SingleFrameEnsemble(InXEnsemble(volume_a)),
            OutXEnsemble(volume_a | volume_b) & LeaveXEnsemble(volume_x),
            SingleFrameEnsemble(InXEnsemble(volume_a | volume_b))
        ])
        return ens


    @staticmethod
    def TISEnsembleSet(volume_a, volume_b, volumes_x, lazy=True):
        myset = []
        for vol in volumes_x:
            myset.append(
                EnsembleFactory.TISEnsemble(volume_a, volume_b, vol, lazy)
            )
        return myset

