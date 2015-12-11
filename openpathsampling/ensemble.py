'''
Created on 03.09.2014

@author: Jan-Hendrik Prinz, David W.H. Swenson
'''

import logging

from openpathsampling.netcdfplus import StorableNamedObject
import openpathsampling as paths

import abc

logger = logging.getLogger(__name__)
init_log = logging.getLogger('openpathsampling.initialization')

# TODO: Make Full and Empty be Singletons to avoid storing them several times!

def join_ensembles(ensemble_list):
    ensemble = None
    for ens in ensemble_list:
        if ensemble is None:
            ensemble = ens
        else:
            ensemble = ensemble | ens
    return ensemble


# note: the cache is not storable, because that would just be silly!
class EnsembleCache(object):
    """Object used by ensembles to enable fast algorithms for basic functions.

    The contents stored in the `can_append`, `can_prepend`, `call`, and
    `check_reverse` dictionaries will depend on the ensemble. Only two of
    these dictionaries should be non-`None` at any time: either the pair
    `call` and `can_append`, or the pair `check_reverse` and `can_prepend`.

    This object also contains basic functions to manage the cache.

    Attributes
    ----------
        start_frame : Snapshot
        prev_last_frame : Snapshot
        direction : +1 or -1
        contents : dictionary
    """
    def __init__(self, direction=None):
        self.start_frame = None
        self.prev_last_frame = None
        self.last_length = None
        self.direction = direction
        self.contents = {}

    def bad_direction_error(self):
        raise RuntimeError("EnsembleCache.direction = " +
                           str(self.direction) + " invalid.") #nocover

    # def clear(self):
    #     self.start_frame = None
    #     self.prev_last_frame = None
    #     self.last_length = None
    #     self.contents = {}

    def check(self, trajectory=None, reset=None):
        """Checks and resets (if necessary) the ensemble cache.
        """
        logger.debug("Checking cache....")
        #logger.debug("traj " + str([id(s) for s in trajectory]))
        logger.debug("start_frame " + str(id(self.start_frame)))
        logger.debug("prev_last " + str(id(self.prev_last_frame)))

        if trajectory is not None:
            # if the first frame has changed, we should reset
            if reset is None:
                lentraj = len(trajectory)
                if self.direction > 0:
                    if trajectory.get_as_proxy(0) != self.start_frame:
                        reset = True
                    else:
                        if lentraj == self.last_length:
                            reset = (trajectory.get_as_proxy(-1) != self.prev_last_frame)
                        elif lentraj == self.last_length + 1:
                            reset = (trajectory.get_as_proxy(-2) != self.prev_last_frame)
                        else:
                            reset = True
                elif self.direction < 0:
                    if trajectory.get_as_proxy(-1) != self.start_frame:
                        reset = True
                    else:
                        if lentraj == self.last_length:
                            reset = (trajectory.get_as_proxy(0) != self.prev_last_frame)
                        elif lentraj == self.last_length + 1:
                            reset = (trajectory.get_as_proxy(1) != self.prev_last_frame)
                        else:
                            reset = True
                else:
                    self.bad_direction_error()
        else:
            reset = True

        self.trusted = not reset 
        self.last_length = len(trajectory)
        if reset:
            logger.debug("Resetting cache " + str(self))
            if self.direction > 0:
                self.start_frame = trajectory.get_as_proxy(0)
                self.prev_last_frame = trajectory.get_as_proxy(-1)
                self.last_length = len(trajectory)
                self.contents = { }
            elif self.direction < 0:
                self.start_frame = trajectory.get_as_proxy(-1)
                self.prev_last_frame = trajectory.get_as_proxy(0)
                self.last_length = len(trajectory)
                self.contents = { }
            else:
                self.bad_direction_error()
        else:
            self.trusted = True
        # by returning reset, we allow the functions that call this to reset
        # other things as well
        if self.direction > 0:
            self.prev_last_frame = trajectory.get_as_proxy(-1)
        elif self.direction < 0:
            self.prev_last_frame = trajectory.get_as_proxy(0)
        else:
            self.bad_direction_error()

        return reset

class Ensemble(StorableNamedObject):
    '''
    Path ensemble object.

    An Ensemble represents a path ensemble, effectively a set of trajectories.
    Typical set operations are allowed, here: and, or, xor, -(without), ~
    (inverse = all - x)     
    
    Examples
    --------    
    >>> EnsembleFactory.TISEnsemble(
    >>>     CVRangeVolume(collectivevariable_A, 0.0, 0.02),
    >>>     CVRangeVolume(collectivevariable_A, 0.0, 0.02),
    >>>     CVRangeVolume(collectivevariable_A, 0.0, 0.08),
    >>>     True
    >>>     )

    Notes
    -----
    Maybe replace - by / to get better notation. So far it has not been used
    '''

    __metaclass__ = abc.ABCMeta

    use_shortcircuit = True

    def __init__(self):
        '''
        A path volume defines a set of paths.
        '''
        super(Ensemble, self).__init__()

    def __eq__(self, other):
        if self is other:
            return True
        return str(self) == str(other)

    @abc.abstractmethod
    def __call__(self, trajectory, trusted=None):
        '''
        Return `True` if the trajectory is part of the path ensemble.

        Parameters
        ----------
        trusted : boolean
            If trusted is not None it overrides the default setting in the ensemble
        '''
        return False

    def check_reverse(self, trajectory, trusted=False):
        return self(trajectory, trusted=False)

    def check(self, trajectory):
        return self(trajectory, trusted = False)

    def trajectory_summary(self, trajectory):
        """
        Return dict with info on how this ensemble "sees" the trajectory.
        """
        return { }

    def trajectory_summary_str(self, trajectory):
        """
        Returns a string with the results of the trajectory_summary function.
        """
        summ = self.trajectory_summary(trajectory)
        if summ == { }:
            return "No summary available"
        else:
            return str(summ)


    def oom_matrix(self, oom):
        """
        Return the oom representation where the OOM is based on a set of volumes

        """

        # Needs to be implemented by the actual class

        return None
    
    def can_append(self, trajectory, trusted=False):
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
        PartInXEnsemble or PartOutXEnsembles. In theory these can only be
        checked if the full range of frames has been generated. This could
        be triggered, when the last frame is reached.  This is even more
        difficult if this depends on the length.
        '''
        return True        

    def can_prepend(self, trajectory, trusted=False):
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
        PartInXEnsemble or PartOutXEnsembles. In theory these can only be checked
        if the full range of frames has been generated. This could be
        triggered, when the last frame is reached.  This is even more
        difficult if this depends on the length.
        '''
        return True        



    def find_valid_slices(self, trajectory, lazy=True, 
                          max_length=None, min_length=1, overlap=1):
        '''
        Return slices (subtrajectories) matching the given ensemble.

        Parameters
        ----------
        trajectory : Trajectory
            the actual trajectory to be splitted into ensemble parts
        lazy : boolean, optional
            if True will use a faster almost linear algorithm, while False
            will run through all possibilities starting with the largest
            ones
        max_length : int > 0, optional
            if set this determines the maximal size to be tested (is mainly
            used in the recursion)
        min_length : int > 0, optional
            if set this determines the minimal size to be tested (in lazy
            mode might no
        overlap : int >= 0, optional
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
                if self.can_append(tt) and end < length:
                    end += 1
                else:
                    if self(tt, trusted=False):
                        ensemble_list.append(slice(start,end))
                        pad = min(overlap, end - start - 1)
                        start = end - pad
                        if end == length:
                            # This means we have reached the end and should stop
                            # All other possible subtraj can only be contained
                            # in already existing ones
                            start = length
                    elif self(tt[0:len(tt)-1], trusted=False):
                        ensemble_list.append(slice(start,end-1))
                        pad = min(overlap, end - start - 2)
                        start = end - pad
                    else:
                        start += 1
                    end = start + min_length

            return ensemble_list

    def split(self, trajectory, lazy=True, max_length=None, min_length=1, overlap=1):
        '''Return list of subtrajectories satisfying the given ensemble.

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

        try:
            # Note here that we use trajectory.lazy() this has the following reason
            # If we would pass the trajectory object itself, then in iterations over
            # snapshots the `for snap in trajectory` will load explicitly the
            # snapshots from storage and so snap is a real Snapshot object.
            # By real I mean that type(snap) is paths.Snapshot equal True!
            # Internally the trajectory just keeps reference objects which are
            # extremely fast to load and since we want the decision to load
            # a snapshots for computing a CV not do always but only if
            # the CV caching decides to we pass trajectory.lazy().
            # The result is that we pass a list of snapshot.proxies for an
            # already stored trajectory and a list of real snapshots for
            # a just created one. The has no speed effect on non-stored
            # trajectories, but makes it faster if the trajectory was loaded or saved
            # and the CV is cached.
            # One more comment, since the idea cannot be always used. The only place
            # where this can fail is if the underlying code uses type(snap) at
            # some point. In this case you need to be able to treat LoaderProxy
            # objects correctly. Since split does not care about the actual snapshots
            # we are safe to use this trick to speed up the evaluation.
            # One last comment about the Proxies. These proxies still behave almost
            # like the real object. If you access any attribute it will be loaded
            # and the actual attribute will be returned. Only difference is operator
            # overloading (which is not used for Snapshots) and type()
            indices = self.find_valid_slices(trajectory.as_proxies(), lazy, max_length,
                                             min_length, overlap)

            return [paths.Trajectory(trajectory[part]) for part in indices]
        except AttributeError:
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
            return UnionEnsemble(self, other)

    def __xor__(self, other):
        if self is other:
            return EmptyEnsemble()
        elif type(other) is EmptyEnsemble:
            return self
        elif type(other) is FullEnsemble:
            return NegatedEnsemble(self)        
        else:
            return SymmetricDifferenceEnsemble(self, other)

    def __and__(self, other):
        if self is other:
            return self
        elif type(other) is EmptyEnsemble:
            return other
        elif type(other) is FullEnsemble:
            return self
        else:
            return IntersectionEnsemble(self, other)

    def __sub__(self, other):
        if self is other:
            return EmptyEnsemble()
        elif type(other) is EmptyEnsemble:
            return self
        elif type(other) is FullEnsemble:
            return EmptyEnsemble()
        else:
            return RelativeComplementEnsemble(self, other)
        
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


class EmptyEnsemble(Ensemble):
    '''
    The empty path ensemble of no trajectories.
    '''
    def __init__(self):
        super(EmptyEnsemble, self).__init__()

    def __call__(self, trajectory, trusted=None):
        return False

    def can_append(self, trajectory, trusted=False):
        return False

    def can_prepend(self, trajectory, trusted=False):
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


class FullEnsemble(Ensemble):
    '''
    The full path ensemble of all possible trajectories.
    '''
    def __init__(self):
        super(FullEnsemble, self).__init__()

    def __call__(self, trajectory, trusted=None):
        return True
    
    def can_append(self, trajectory, trusted=False):
        return True

    def can_prepend(self, trajectory, trusted=False):
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


class NegatedEnsemble(Ensemble):
    '''
    Negates an Ensemble and simulates a `not` statement
    '''
    def __init__(self, ensemble):
        super(NegatedEnsemble, self).__init__()
        self.ensemble = ensemble
        
    def __call__(self, trajectory, trusted=None):
        return not self.ensemble(trajectory, trusted)

    def can_append(self, trajectory, trusted=False):
        # We cannot guess the result here so keep on running forever
        return True

    def can_prepend(self, trajectory, trusted=False):
        # We cannot guess the result here so keep on running forever
        return True

    def __str__(self):
        return 'not ' + str(self.ensemble)


class EnsembleCombination(Ensemble):
    '''
    Logical combination of two ensembles
    '''
    def __init__(self, ensemble1, ensemble2, fnc, str_fnc):
        super(EnsembleCombination, self).__init__()
        self.ensemble1 = ensemble1
        self.ensemble2 = ensemble2
        self.fnc = fnc
        self.sfnc = str_fnc

    def to_dict(self):
        return { 'ensemble1' : self.ensemble1, 'ensemble2' : self.ensemble2 }

    def __call__(self, trajectory, trusted=None):
        # Shortcircuit will automatically skip the second part of the combination if the result does not depend on it!
        # This makes sense since the expensive part is the ensemble testing not computing two logic operations
        if Ensemble.use_shortcircuit:
            a = self.ensemble1(trajectory, trusted)
            logger.debug("Combination is " + self.__class__.__name__)
            logger.debug("Combination: " + self.ensemble1.__class__.__name__ + 
                         " is "+str(a))
            logger.debug("Combination: " + self.ensemble2.__class__.__name__ +
                         " is " +str(self.ensemble2(trajectory, trusted)))
            logger.debug("Combination: returning " + 
                         str(self.fnc(a,self.ensemble2(trajectory,trusted))))
            res_true = self.fnc(a, True)
            res_false = self.fnc(a, False)
            if res_false == res_true:
                # result is independent of ensemble_b so ignore it
                return res_true
            else:
                b = self.ensemble2(trajectory, trusted)
                return self.fnc(a, b)
        else:
            return self.fnc(self.ensemble1(trajectory, trusted), self.ensemble2(trajectory, trusted))

    # Forward / Backward is tricky
    # We can do the following. If a or b is true this means that the real
    # result could be false or true, we just keep going but we should have
    # stopped. If a or b is false this means for that ensemble continuing is
    # not feasible and so false really means false. To check if a logical
    # combination should be continued just try for all true values a
    # potential false and check if we should continue.

    def _continue_fnc(self, a, b):
        fnc = self.fnc
        res = fnc(a,b)
        if a == True:
            res |= fnc(False, b)
        if b == True:
            res |= fnc(a, False)
        if a == True and b == True:
            res |= fnc(False, False)

        return res

    def can_append(self, trajectory, trusted=False):
        if Ensemble.use_shortcircuit:
            a = self.ensemble1.can_append(trajectory, trusted)
            logger.debug("Combination is " + self.__class__.__name__)
            logger.debug("Combination.can_append: " + 
                         self.ensemble1.__class__.__name__ + " is "+str(a))
            logger.debug("Combination.can_append: " 
                         + self.ensemble2.__class__.__name__ +
                         " is " +
                         str(self.ensemble2.can_append(trajectory, trusted)))
            logger.debug("Combination.can_append: returning " + 
                         str(self.fnc(a,self.ensemble2.can_append(trajectory,trusted))))
            res_true = self._continue_fnc(a, True)
            res_false = self._continue_fnc(a, False)
            if res_false == res_true:
                # result is independent of ensemble_b so ignore it
                return res_true
            else:
                b = self.ensemble2.can_append(trajectory, trusted)
                #logger.debug("b is " + str(b))
                if b == True:
                    #logger.debug("Will return res_true")
                    return res_true
                else:
                    #logger.debug("Will return res_false")
                    return res_false
        else:
            return self.fnc(self.ensemble1.can_append(trajectory, trusted), 
                            self.ensemble2.can_append(trajectory, trusted))

    def can_prepend(self, trajectory, trusted=False):
        if Ensemble.use_shortcircuit:
            a = self.ensemble1.can_prepend(trajectory, trusted)
            res_true = self._continue_fnc(a, True)
            res_false = self._continue_fnc(a, False)
            if res_false == res_true:
                # result is independent of ensemble_b so ignore it
                return res_true
            else:
                b = self.ensemble2.can_prepend(trajectory, trusted)
                if b == True:
                    return res_true
                else:
                    return res_false
        else:
            return self.fnc(self.ensemble1.can_prepend(trajectory, trusted), 
                            self.ensemble2.can_prepend(trajectory, trusted))

    def __str__(self):
#        print self.sfnc, self.ensemble1, self.ensemble2, self.sfnc.format('(' + str(self.ensemble1) + ')' , '(' + str(self.ensemble1) + ')')
        return self.sfnc.format('(\n' + Ensemble._indent(str(self.ensemble1)) + '\n)' , '(\n' + Ensemble._indent(str(self.ensemble2)) + '\n)')


class UnionEnsemble(EnsembleCombination):
    def __init__(self, ensemble1, ensemble2):
        super(UnionEnsemble, self).__init__(ensemble1, ensemble2, fnc = lambda a,b : a or b, str_fnc = '{0}\nor\n{1}')


class IntersectionEnsemble(EnsembleCombination):
    def __init__(self, ensemble1, ensemble2):
        super(IntersectionEnsemble, self).__init__(ensemble1, ensemble2, fnc = lambda a,b : a and b, str_fnc = '{0}\nand\n{1}')


class SymmetricDifferenceEnsemble(EnsembleCombination):
    def __init__(self, ensemble1, ensemble2):
        super(SymmetricDifferenceEnsemble, self).__init__(ensemble1, ensemble2, fnc = lambda a,b : a ^ b, str_fnc = '{0}\nxor\n{1}')


class RelativeComplementEnsemble(EnsembleCombination):
    def __init__(self, ensemble1, ensemble2):
        super(RelativeComplementEnsemble, self).__init__(ensemble1, ensemble2, fnc = lambda a,b : a and not b, str_fnc = '{0}\nand not\n{1}')


class SequentialEnsemble(Ensemble):
    """Ensemble which satisfies several subensembles in sequence.

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

        self._use_cache = True # cache can be turned off
        self._cache_can_append = EnsembleCache(+1)
        self._cache_call = EnsembleCache(+1)
        self._cache_can_prepend = EnsembleCache(-1)
        self._cache_check_reverse = EnsembleCache(-1)

        # sanity checks
        if len(self.min_overlap) != len(self.max_overlap):
            raise ValueError("len(min_overlap) != len(max_overlap)")
        if len(self.min_overlap) != len(self.ensembles)-1:
            raise ValueError("Number of overlaps doesn't match number of transitions")
        for i in range(len(self.min_overlap)):
            if min_overlap[i] > max_overlap[i]:
                raise ValueError("min_overlap greater than max_overlap!")

    def update_cache(self, cache, ens_num, ens_from, subtraj_from):
        """Updates the given cache.

        Parameters
        ----------
        ens_num : integer
            current value of `ens_num` in the sequential ensemble
        ens_from : integer
            current "start" ensemble index. For forward-direction caches,
            this is ens_first. For reverse-direction caches, this is
            ens_final. The "initial" (in the appropriate direction) frame is
            assigned to this ensemble
        subtraj_from : integer
            index of the "start" frame of the subtrajectory in this
            subensemble. For forward-direction caches, this is the first
            frame of the subtrajectory. For reverse-direction caches, this
            is the final frame of the subtrajectory.
        """
        if ens_num == "keep":
            ens_num = cache.contents['ens_num']
        if ens_from == "keep":
            ens_from = cache.contents['ens_from']
        if subtraj_from == "keep":
            subtraj_from = cache.contents['subtraj_from']

        cache.contents['ens_num'] = ens_num
        cache.contents['ens_from'] = ens_from
        cache.contents['subtraj_from'] = subtraj_from
        logger.debug("Setting cache | ens_num " + str(ens_num) +
                     " | ens_from " + str(ens_from) +
                     " | subtraj_from " + str(subtraj_from)
                    )
        logger.debug("Cache is Trusted: " + str(cache.trusted))

    def assign_frames(self, cache, ens_num, subtraj_first=None, subtraj_final=None):
        if ens_num is None:
            cache.contents['assignments'] = { }
        else:
            cache.contents['assignments'][ens_num] = slice(subtraj_first, subtraj_final) 
        logger.debug("Cache assignments: " + str(cache.contents['assignments']))


    def transition_frames(self, trajectory, trusted=None):
        # it is easiest to understand this decision tree as a simplified
        # version of the can_append decision tree; see that for detailed
        # comments
        #self._check_cache(trajectory, function="call")

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


    def __call__(self, trajectory, trusted=None):
        logger.debug("Looking for transitions in trajectory " + str(trajectory))
        transitions = self.transition_frames(trajectory, trusted)
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

    def _find_subtraj_final(self, traj, subtraj_first, ens_num,
                            last_checked=None):
        """
        Find the longest subtrajectory of trajectory which starts at
        subtraj_first and satifies self.ensembles[ens_num].can_append

        Returns
        -------
        int
            Frame of traj which is the final frame for a subtraj starting at
            subtraj_first and satisfying self.ensembles.can_append[ens_num]
        """
        if last_checked is None:
            subtraj_final = subtraj_first
        else:
            subtraj_final = max(last_checked, subtraj_first)
        traj_final = len(traj)
        ens = self.ensembles[ens_num]
        subtraj = traj[slice(subtraj_first, subtraj_final+1)]
        # if we're in the ensemble or could eventually be in the ensemble,
        # we keep building the subtrajectory

        # TODO: this doesn't actually reflect the cleanest behavior: should
        # be the proper hybrid definition where we can append until/unless
        # we overshoot
        logger.debug("*Traj slice " + str(subtraj_first) + " " + 
                     str(subtraj_final+1) + " / " + str(traj_final))
        #logger.debug("Ensemble " + str(ens.__class__.__name__))# + str(ens))
        #logger.debug("Can-app " + str(ens.can_append(subtraj, trusted=True)))
        #logger.debug("Call    " + str(ens(subtraj, trusted=True)))
        while ( (ens.can_append(subtraj, trusted=True) or 
                 ens(subtraj, trusted=True)
                ) and subtraj_final < traj_final):
            subtraj_final += 1
            # TODO: replace with append; probably faster
            subtraj = traj[slice(subtraj_first, subtraj_final+1)]
            logger.debug(" Traj slice " + str(subtraj_first) + " " + 
                         str(subtraj_final+1) + " / " + str(traj_final))
        return subtraj_final
    
    def _find_subtraj_first(self, traj, subtraj_final, ens_num,
                            last_checked=None):
        if last_checked is None:
            subtraj_first = subtraj_final-1
        else:
            subtraj_first = min(last_checked, subtraj_final-1)
        traj_first = 0
        ens = self.ensembles[ens_num]
        subtraj = traj[slice(subtraj_first, subtraj_final)]
        logger.debug("*Traj slice " + str(subtraj_first) + " " + 
                     str(subtraj_final) + " / " + str(len(traj)))
        #logger.debug("Ensemble " + str(ens.__class__.__name__))# + str(ens))
        #logger.debug("Can-app " + str(ens.can_prepend(subtraj, trusted=True)))
        #logger.debug("Call    " + str(ens(subtraj, trusted=True)))
        while ( (ens.can_prepend(subtraj, trusted=True) or 
                 ens.check_reverse(subtraj, trusted=True)
                ) and subtraj_first >= traj_first):
            subtraj_first -= 1
            subtraj = traj[slice(subtraj_first, subtraj_final)]
            logger.debug(" Traj slice " + str(subtraj_first+1) + " " + 
                         str(subtraj_final) + " / " + str(len(traj)))
        return subtraj_first+1


    def can_append(self, trajectory, trusted=False):
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
        #self._check_cache(trajectory, function="can_append")
        cache = self._cache_can_append
        if trusted:
            cache.trusted = True

        subtraj_first = 0
        ens_num = 0
        ens_first = 0


        if self._use_cache: 
            cache.check(trajectory)
            if cache.contents == { }:
                self.update_cache(cache, 0, 0, 0)
                self.assign_frames(cache, None, None, None)
            else:
                subtraj_first = cache.contents['subtraj_from']
                ens_num = cache.contents['ens_num']
                ens_first = cache.contents['ens_from']


        traj_final = len(trajectory)
        final_ens = len(self.ensembles)-1
        #print traj_final, final_ens
        # logging startup
        logger.debug(
            "Beginning can_append with subtraj_first=" + str(subtraj_first)
            + "; ens_first=" + str(ens_first) + "; ens_num=" + str(ens_num)
        )
        logger.debug(
            "Can-append sees a trusted cache: " + str(cache.trusted)
        )
        for i in range(len(self.ensembles)):
            ens = self.ensembles[i]
            logger.debug("Ensemble " + str(i) + " : " + ens.__class__.__name__)

        while True: #  main loop, with various 
            if self._use_cache and cache.trusted:
                last_checked = trajectory.index(cache.prev_last_frame)-1
            else:
                last_checked = None
            subtraj_final = self._find_subtraj_final(
                trajectory, subtraj_first, ens_num, last_checked
            )
            logger.debug(
                "Subtraj for ens " + str(ens_num) + " : " +
                "("+str(subtraj_first)+","+str(subtraj_final)+")"
            )
            if subtraj_final - subtraj_first > 0:
                subtraj = trajectory[slice(subtraj_first, subtraj_final)]
                if ens_num == final_ens:
                    if subtraj_final == traj_final:
                        # we're in the last ensemble and the whole
                        # trajectory is assigned: can we append?
                        ens = self.ensembles[ens_num]
                        logger.debug("Returning can_append for " 
                                     + str(ens.__class__.__name__))
                        self.update_cache(cache, ens_num, 
                                          ens_first, subtraj_first)
                        return ens.can_append(subtraj, trusted=True)
                    else:
                        logger.debug(
                            "Returning false due to incomplete assigns: " + 
                            str(subtraj_final) + "!=" + str(traj_final)
                        )
                        return False # in final ensemble, not all assigned
                else:
                    # subtraj existed, but not yet final ensemble
                    # so we start with the next ensemble
                    if subtraj_final != traj_final and not self.ensembles[ens_num](subtraj, trusted=cache.trusted):
                        logger.debug(
                            "Couldn't assign frames " + str(subtraj_first) +
                            " through " + str(subtraj_final) + 
                            " to ensemble " + str(ens_num) + ": No match"
                        )
                    else:
                        logger.debug(
                            "Assigning frames " + str(subtraj_first) +
                            " through " + str(subtraj_final) + 
                            " to ensemble " + str(ens_num)
                        )
                        self.assign_frames(cache, ens_num, subtraj_first,
                                           subtraj_final)
                        self.update_cache(cache, ens_num, ens_first, subtraj_first)
                    ens_num += 1
                    subtraj_first = subtraj_final
                    logger.debug("Moving to the next ensemble " + str(ens_num))
            else: 
                if subtraj_final == traj_final:
                    # all frames assigned, but not all ensembles finished;
                    # next frame might satisfy next ensemble
                    if self._use_cache:
                        prev_slice = cache.contents['assignments'][ens_num-1]
                        prev_subtraj = trajectory[prev_slice]
                        prev_ens = self.ensembles[ens_num-1]
                        if prev_ens.can_append(prev_subtraj, trusted=True):
                            logger.debug(
                                "Premature promotion: returning to ensemble " + 
                                str(ens_num-1)
                            )
                            ens_num -= 1
                            subtraj_first = "keep"
                        
                        self.update_cache(cache, ens_num, ens_first, subtraj_first)
                    logger.debug("All frames assigned, more ensembles to go: returning True")
                    return True
                elif self.ensembles[ens_num](paths.Trajectory([])):
                    logger.debug("Moving on because of allowed zero-length ensemble")
                    ens_num += 1
                    subtraj_first = subtraj_final
                    self.update_cache(cache, ens_num, ens_first, subtraj_first)
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
                        self.update_cache(cache, ens_num, ens_first, subtraj_first)


    def can_prepend(self, trajectory, trusted=False):
        # based on .can_append(); see notes there for algorithm details
        cache = self._cache_can_prepend
        if trusted:
            cache.trusted = True

        traj_first = 0
        first_ens = 0
        subtraj_final = len(trajectory)
        ens_final = len(self.ensembles)-1
        ens_num = ens_final

        if self._use_cache:
            cache.check(trajectory)
            if cache.contents == { }:
                self.update_cache(cache, ens_num, first_ens, subtraj_final)
                self.assign_frames(cache, None, None, None)
            else:
                logger.debug("len(traj)="+str(len(trajectory)) + 
                             "cache_from="+str(cache.contents['subtraj_from']))
                subtraj_from = cache.contents['subtraj_from']
                if subtraj_from == None:
                    subtraj_from = 0
                subtraj_final = len(trajectory)+subtraj_from
                ens_num = cache.contents['ens_num']
                ens_final = cache.contents['ens_from']


        # logging startup
        logger.debug("Beginning can_prepend with ens_num:" + str(ens_num) +
                     "  ens_final:" + str(ens_final) + "  subtraj_final " +
                     str(subtraj_final)
                    )
        for i in range(len(self.ensembles)):
            logger.debug(
                "Ensemble " + str(i) + 
                " : " + self.ensembles[i].__class__.__name__
            )

        while True:
            if self._use_cache and cache.trusted:
                last_checked = trajectory.index(cache.prev_last_frame)+1
            else:
                last_checked = None
            subtraj_first = self._find_subtraj_first(
                trajectory, subtraj_final, ens_num, last_checked)

            assign_final = subtraj_final - len(trajectory)
            if assign_final == 0:
                assign_final = None
            logger.debug(
                str(ens_num) + " : " +
                "("+str(subtraj_first)+","+str(subtraj_final)+")"
            )
            if subtraj_final - subtraj_first > 0:
                subtraj = trajectory[slice(subtraj_first, subtraj_final)]
                if ens_num == first_ens:
                    if subtraj_first == traj_first:
                        logger.debug("Returning can_prepend")
                        self.update_cache(cache, ens_num, ens_final,
                                          assign_final)
                        return self.ensembles[ens_num].can_prepend(subtraj,
                                                                  trusted=True)
                    else:
                        logger.debug(
                            "Returning false due to incomplete assigns: " + 
                            str(subtraj_first) + "!=" + str(traj_first)
                        )
                        return False
                else:
                    if subtraj_first != traj_first and not self.ensembles[ens_num](subtraj, trusted=True):
                        logger.debug(
                            "Couldn't assign frames " + str(subtraj_first) +
                            " through " + str(subtraj_final) + 
                            " to ensemble " + str(ens_num) + ": No match"
                        )
                    else:
                        logger.debug(
                            "Assigning frames " + str(subtraj_first) +
                            " through " + str(subtraj_final) + 
                            " to ensemble " + str(ens_num)
                        )
                        assign_first = subtraj_first-len(trajectory)
                        self.assign_frames(cache, ens_num, assign_first,
                                           assign_final)
                        self.update_cache(cache, ens_num, ens_final,
                                          assign_final)
                    ens_num -= 1
                    subtraj_final = subtraj_first
                    logger.debug("Moving to the next ensemble " + str(ens_num))
            else:
                if subtraj_first == traj_first:
                    if self._use_cache:
                        prev_slice = cache.contents['assignments'][ens_num+1]
                        logger.debug("prev_slice " + str(prev_slice))
                        prev_subtraj = trajectory[prev_slice]
                        logger.debug("prev_subtraj " + str(prev_subtraj))
                        logger.debug("traj " + str(trajectory))
                        prev_ens = self.ensembles[ens_num+1]
                        if prev_ens.can_prepend(prev_subtraj, trusted=True):
                            logger.debug(
                                "Premature promotion: returning to ensemble " + 
                                str(ens_num+1)
                            )
                            ens_num += 1
                            assign_final = "keep"

                        logger.debug("(first, final)" + str( (subtraj_first,
                                                              subtraj_final)))
                        self.update_cache(cache, ens_num, ens_final, 
                                          assign_final)
                    logger.debug("All frames assigned, more ensembles to go: returning True")
                    return True
                elif self.ensembles[ens_num](paths.Trajectory([])):
                    logger.debug("Moving on because of allowed zero-length ensemble")
                    ens_num -= 1
                    subtraj_final = subtraj_first
                    self.update_cache(cache, ens_num, ens_final, subtraj_final)
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
                        self.update_cache(cache, ens_num, ens_final, subtraj_final)

    def __str__(self):
        head = "[\n"
        tail = "\n]"
        sequence_str = ",\n".join([str(ens) for ens in self.ensembles])
        return head+sequence_str+tail



class LengthEnsemble(Ensemble):
    '''
    The ensemble of trajectories of a given length
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
    
    def __call__(self, trajectory, trusted=None):
        length = len(trajectory)
        if type(self.length) is int:
            return length == self.length
        else:
            return length >= self.length.start and (self.length.stop is None or length < self.length.stop)
        
    def can_append(self, trajectory, trusted=False):
        length = len(trajectory)
        if type(self.length) is int:
            return_value = (length < self.length)
            logger.debug("LengthEnsemble.can_append: Segment length " + 
                         str(length) + " < " + str(self.length) + " : " +
                         str(return_value))
            return return_value
        else:
            return self.length.stop is None or length < self.length.stop - 1

    def can_prepend(self, trajectory, trusted=False):
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


class VolumeEnsemble(Ensemble):
    '''
    Path ensembles based on the Volume object
    '''    
    def __init__(self, volume, trusted = True):
        super(VolumeEnsemble, self).__init__()
        self.volume = volume
        self.trusted = trusted

    @property
    def _volume(self):
        '''
        The volume that is used in the specification.
        '''
        return self.volume


class AllInXEnsemble(VolumeEnsemble):
    '''
    Ensemble of trajectories with all frames in the given volume
    '''

    def can_append(self, trajectory, trusted=False):
        if len(trajectory) == 0:
            return True
        if trusted == True:
            return self(trajectory[slice(len(trajectory)-1, None)], trusted)
        else:
            return self(trajectory)

    def can_prepend(self, trajectory, trusted=False):
        if len(trajectory) == 0:
            return True
        if trusted == True:
            return self(trajectory[slice(0,1)], trusted)
        else:
            return self(trajectory)
        
    
    def __call__(self, trajectory, trusted=None):
        if len(trajectory) == 0:
            return False
        if trusted == True:
            #print "trusted"
            frame = trajectory.get_as_proxy(-1)
            return self._volume(frame)
        else:
            #logger.debug("Calling volume untrusted "+repr(self))
            for frame in trajectory.as_proxies():
                if not self._volume(frame):
                    return False
            return True

    def check_reverse(self, trajectory, trusted=False):
        # order in this one only matters if it is trusted
        if trusted:
            #print "Rev Trusted"
            frame = trajectory.get_as_proxy(0)
            return self._volume(frame)
        else:
            #print "Rev UnTrusted"
            return self(trajectory) # in this case, order wouldn't matter


    def __invert__(self):
        return PartOutXEnsemble(self.volume, self.frames, self.trusted)

    def __str__(self):
        return 'x[t] in {0} for all t'.format(self._volume)



class AllOutXEnsemble(AllInXEnsemble):
    '''
    Ensemble of trajectories with all frames outside the given volume
    '''    
    @property
    def _volume(self):
        return ~ self.volume
    
    def __str__(self):
        return 'x[t] in {0} for all t'.format(self._volume)

    def __invert__(self):
        return PartInXEnsemble(self.volume, self.frames, self.trusted)


class PartInXEnsemble(VolumeEnsemble):
    '''
    Ensemble of trajectory with at least one frame in the volume
    '''

    def __str__(self):
        return 'exists t such that x[t] in {0}'.format(self._volume)

    def __call__(self, trajectory, trusted=None):
        '''
        Returns True if the trajectory is part of the PathEnsemble
        
        Parameters
        ----------
        trajectory : Trajectory
            The trajectory to be checked
        '''
        for frame in trajectory.as_proxies():
            if self._volume(frame):
                return True
        return False

    def __invert__(self):
        return AllOutXEnsemble(self.volume, self.frames, self.trusted)


class PartOutXEnsemble(PartInXEnsemble):
    '''
    Ensemble of trajectories with at least one frame outside the volume
    '''
    def __str__(self):
        return 'exists t such that x[t] in {0}'.format(self._volume)
      
    @property
    def _volume(self):
        # effectively use PartInXEnsemble but with inverted volume
        return ~ self.volume

    def __invert__(self):
        return AllInXEnsemble(self.volume, self.frames, self.trusted)

    def __call__(self, trajectory, trusted=None):
        for frame in trajectory.as_proxies():
            if self._volume(frame):
                return True
        return False


class ExitsXEnsemble(VolumeEnsemble):
    """
    Represents an ensemble where two successive frames from the selected
    frames of the trajectory crossing from inside to outside the given volume.
    """
    def __init__(self, volume, trusted=False):
        # changing the defaults for frames and trusted; prevent single frame
        super(ExitsXEnsemble, self).__init__(volume, trusted)

    def __str__(self):
        domain = 'exists x[t], x[t+1] '
        result = 'such that x[t] in {0} and x[t+1] not in {0}'.format(
                            self._volume)
        return domain+result

    def __call__(self, trajectory, trusted=None):
        subtraj = trajectory
        for i in range(len(subtraj)-1):
            frame_i = subtraj.get_as_proxy(i)
            if self._volume(frame_i):
                frame_iplus = subtraj.get_as_proxy(i+1)
                if not self._volume(frame_iplus):
                    return True
        return False


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

    def __call__(self, trajectory, trusted=None):
        subtraj = trajectory
        for i in range(len(subtraj)-1):
            frame_i = subtraj.get_as_proxy(i)
            if not self._volume(frame_i):
                frame_iplus = subtraj.get_as_proxy(i+1)
                if self._volume(frame_iplus):
                    return True
        return False


class WrappedEnsemble(Ensemble):
    '''
    Wraps an ensemble to alter it or the way it sees a trajectory
    '''
    def __init__(self, ensemble):
        super(WrappedEnsemble, self).__init__()
        self.ensemble = ensemble

        # you can also build wrapped ensembles with more flexibility when using
        # a property for _new_ensemble
        self._new_ensemble = self.ensemble
        self.trusted = None
        self._cache_can_append = EnsembleCache(+1)
        self._cache_call = EnsembleCache(+1)
        self._cache_can_prepend = EnsembleCache(+1) #TODO: this is weird
        # cache_can_prepend has to think it is going forward because the
        # frames given to it are from a forward growing trajectory... only
        # later is everything turned around

    def __call__(self, trajectory, trusted=None):
        return self._new_ensemble(self._alter(trajectory), trusted)

    def _alter(self, trajectory):
        return trajectory
        
    def can_append(self, trajectory, trusted=None):
        return self._new_ensemble.can_append(self._alter(trajectory),
                                             trusted)

    def can_prepend(self, trajectory, trusted=None):
        return self._new_ensemble.can_prepend(self._alter(trajectory))


class SlicedTrajectoryEnsemble(WrappedEnsemble):
    '''
    Alters trajectories given as arguments by taking Python slices.
    '''
    def __init__(self, ensemble, region):
        super(SlicedTrajectoryEnsemble, self).__init__(ensemble)
        if type(region) == int:
            if region == -1:
                self.region = slice(region,None)
            else:
                self.region = slice(region, region+1)
        else:
            self.region = region

    def _alter(self, trajectory):
        return trajectory[self.region]

    def __str__(self):
        # TODO: someday may add different string support for slices with
        # only one frame
        start = "" if self.region.start is None else str(self.region.start)
        stop = "" if self.region.stop is None else str(self.region.stop)
        step = "" if self.region.step is None else " every "+str(self.region.step)
        return ("(" + self.ensemble.__str__() +
                " in {" + start + ":" + stop + "}" + step + ")")


class SuffixTrajectoryEnsemble(WrappedEnsemble):
    '''
    Ensemble which prepends its trajectory to a given trajectory.

    Used in backward shooting.
    '''
    def __init__(self, ensemble, add_trajectory):
        super(SuffixTrajectoryEnsemble, self).__init__(ensemble)
        self.add_trajectory = add_trajectory
        self._cached_trajectory = paths.Trajectory(add_trajectory.as_proxies())

    def _alter(self, trajectory):
        logger.debug("Starting Suffix._alter")
        #logger.debug("altered " + str([id(i) for i in self._cached_trajectory]))
        reset = self._cache_can_prepend.check(trajectory)
        #logger.debug("altered " + str([id(i) for i in self._cached_trajectory]))
        #logger.debug("traj    " + str([id(i) for i in trajectory]))
        #logger.debug("trajrev " + str([id(i) for i in trajectory.reversed]))
        #reset = False
        if not reset:
            logger.debug("BackwardPrended was not reset")
            first_frame = trajectory.get_as_proxy(-1)
            if self._cached_trajectory.get_as_proxy(0) != first_frame:
                self._cached_trajectory.insert(0, first_frame)
        else:
            self._cached_trajectory = trajectory.reversed + self.add_trajectory

        #logger.debug("revtraj " + str([id(i) for i in revtraj]))
        #logger.debug("add     " + str([id(i) for i in self.add_trajectory]))
        #logger.debug("altered " + str([id(i) for i in self._cached_trajectory]))

        return self._cached_trajectory

    def can_append(self, trajectory, trusted=None):
        raise RuntimeError("SuffixTrajectoryEnsemble.can_append is nonsense.")


class PrefixTrajectoryEnsemble(WrappedEnsemble):
    '''
    Ensemble which appends its trajectory to a given trajectory.

    Used in forward shooting.
    '''
    def __init__(self, ensemble, add_trajectory):
        super(PrefixTrajectoryEnsemble, self).__init__(ensemble)
        self.add_trajectory = add_trajectory
        self._cached_trajectory = paths.Trajectory(add_trajectory.as_proxies())

    def _alter(self, trajectory):
        logger.debug("Starting _alter")
        reset = self._cache_can_append.check(trajectory)
        if not reset:
            final_frame = trajectory.get_as_proxy(-1)
            if self._cached_trajectory.get_as_proxy(-1) != final_frame:
                self._cached_trajectory.append(final_frame)
        else: 
            logger.debug("doing it oldstyle")
            self._cached_trajectory = self.add_trajectory + trajectory

        # DEBUG 
        #logger.debug("add   " + str([i for i in self.add_trajectory]))
        #logger.debug("traj  " + str([i for i in trajectory]))
        #logger.debug("cache " + str([i for i in self._cached_trajectory]))
        #oldstyle = self.add_trajectory + trajectory
        #for (t,b) in zip(self._cached_trajectory,
                         #self.add_trajectory+trajectory):
            #logger.debug(str(t) + " ?=? " + str(b))
            #assert(t == b)
        #assert(len(self._cached_trajectory) == len(oldstyle))

        return self._cached_trajectory

    def can_prepend(self, trajectory, trusted=None):
        raise RuntimeError("PrefixTrajectoryEnsemble.can_prepend is nonsense.")


class ReversedTrajectoryEnsemble(WrappedEnsemble):
    '''
    Ensemble based on reversing the trajectory.
    '''
    def _alter(self, trajectory):
        return trajectory.reverse()


class AppendedNameEnsemble(WrappedEnsemble):
    '''
    Add string to ensemble name: allows multiple copies of an ensemble.
    '''
    def __init__(self, ensemble, label):
        self.label = label
        super(AppendedNameEnsemble, self).__init__(ensemble)

    def __str__(self):
        return self.ensemble.__str__() + " " + self.label



class OptionalEnsemble(WrappedEnsemble):
    '''
    An ensemble which is optional for SequentialEnsembles.
    '''

    def __init__(self, ensemble):
        super(OptionalEnsemble, self).__init__(ensemble)
        self._new_ensemble = LengthEnsemble(0) | self.ensemble

    def __str__(self):
        return "{"+self.ensemble.__str__()+"} (OPTIONAL)"


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


class MinusInterfaceEnsemble(SequentialEnsemble):
    '''
    This creates an ensemble for the minus interface. 

    Parameters
    ----------
    state_vol : Volume
        The Volume which defines the state for this minus interface
    innermost_vols : list of Volume
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
    T.S. van Erp. Phys. Rev. Lett.

    D.W.H. Swenson and P.G. Bolhuis. J. Chem. Phys. 141, 044101 (2014). 
    doi:10.1063/1.4890037
    '''

    # don't store unnecessary stuff we recreate at initialization
    # TODO: Check with David if it makes sense to store these and allow
    # them being used in __init__ instead of the self-made ones

    _excluded_attr = ['ensembles', 'min_overlap', 'max_overlap']

    def __init__(self, state_vol, innermost_vols, n_l=2, greedy=False):
        if (n_l < 2):
            raise ValueError("The number of segments n_l must be at least 2")

        self.state_vol = state_vol
        try:
            innermost_vols = list(innermost_vols)
        except TypeError:
            innermost_vols = [innermost_vols]

        self.innermost_vols = innermost_vols
        #self.innermost_vol = paths.volume.join_volumes(self.innermost_vols)
        self.innermost_vol = paths.FullVolume()
        for vol in self.innermost_vols:
            self.innermost_vol = self.innermost_vol & vol
        self.greedy = greedy
        inA = AllInXEnsemble(state_vol)
        outA = AllOutXEnsemble(state_vol)
        outX = AllOutXEnsemble(self.innermost_vol)
        inX = AllInXEnsemble(self.innermost_vol)
        leaveX = PartOutXEnsemble(self.innermost_vol)
        interstitial = outA & inX
        segment_ensembles = [paths.TISEnsemble(state_vol, state_vol, inner)
                             for inner in self.innermost_vols]

        self._segment_ensemble = join_ensembles(segment_ensembles)

        #interstitial = AllInXEnsemble(self.innermost_vol - state_vol)
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

    def populate_minus_ensemble(self, partial_traj, minus_replica_id, engine):
        """
        Generate a sample for the minus ensemble by extending `partial_traj`

        Parameters
        ----------
        partial_traj : Trajectory
            trajectory to extend
        minus_replica_id : integer or string
            replica ID for this sample
        engine : DynamicsEngine
            engine to use for MD extension
        """
        last_frame = partial_traj[-1]
        if not self._segment_ensemble(partial_traj):
            raise RuntimeError(
                "Invalid input trajectory for minus extension. (Not A-to-A?)"
            )
        extension = engine.generate(last_frame,
                                    [self.can_append])
        first_minus = paths.Trajectory(partial_traj + extension[1:])
        minus_samp = paths.Sample(
            replica=minus_replica_id,
            trajectory=first_minus,
            ensemble=self
        )
        logger.info(first_minus.summarize_by_volumes_str(
            {"A" : self.state_vol,
             "I" : ~self.state_vol & self.innermost_vol,
             "X" : ~self.innermost_vol})
        )
        return minus_samp

class TISEnsemble(SequentialEnsemble):
    """An ensemble for TIS (or AMS).

    Begin in `initial_states`, end in either `initial_states` or
    `final_states`, and cross `interface`.

    Attributes
    ----------
    initial_states : Volume or list of Volume
        Volume(s) that only the first or last frame may be in
    final_states : Volume or list of Volume
        Volume(s) that only the last frame may be in
    interface : Volume
        Volume which the trajectory must exit to be accepted
    orderparameter : CollectiveVariable
        CV to be used as order parameter for this
    """
    def __init__(self, initial_states, final_states, interface,
                 orderparameter=None):
        # regularize to list of volumes
        # without orderparameter, some info can't be obtained
        try:
            n_initial_states = len(initial_states)
        except TypeError:
            n_initial_states = 1
            initial_states = [initial_states]

        try:
            n_final_states = len(final_states)
        except TypeError:
            n_final_states = 1
            final_states = [final_states]

        volume_a = paths.volume.join_volumes(initial_states)
        volume_b = paths.volume.join_volumes(final_states)

        super(TISEnsemble, self).__init__([
            AllInXEnsemble(volume_a) & LengthEnsemble(1),
            AllOutXEnsemble(volume_a | volume_b) & PartOutXEnsemble(interface),
            AllInXEnsemble(volume_a | volume_b) & LengthEnsemble(1)
        ])

        self.initial_states = initial_states
        self.final_states = final_states
        self.interface = interface
#        self.name = interface.name
        self.orderparameter = orderparameter

    def trajectory_summary(self, trajectory):
        initial_state_i = None
        final_state_i = None
        for state_i in range(len(self.initial_states)):
            if self.initial_states[state_i](trajectory.get_as_proxy(0)):
                initial_state_i = state_i
                break
        all_states = self.initial_states + self.final_states
        for state_i in range(len(all_states)):
            if all_states[state_i](trajectory.get_as_proxy(-1)):
                final_state_i = state_i
                break

        if self.orderparameter is not None:
            lambda_traj = self.orderparameter(trajectory)
            min_lambda = min(lambda_traj)
            max_lambda = max(lambda_traj)
        else:
            min_lambda = None
            max_lambda = None

        return {
            'initial_state' : initial_state_i,
            'final_state' : final_state_i,
            'max_lambda' : max_lambda,
            'min_lambda' : min_lambda
        }


    def trajectory_summary_str(self, trajectory):
        summ = self.trajectory_summary(trajectory)
        all_states = self.initial_states + self.final_states
        # TODO: remove the .name from this when string returns correctly
        init_st_i = summ['initial_state']
        fin_st_i = summ['final_state']
        # TODO: how can we have None?
        if init_st_i == None:
            init_st = "None"
        else:
            init_st = str(self.initial_states[summ['initial_state']].name)
        if fin_st_i == None:
            fin_st = "None"
        else:
            fin_st = str(all_states[summ['final_state']].name)

        if self.orderparameter is not None:
            opname = self.orderparameter.name
        else:
            opname = "None"
        min_l = str(summ['min_lambda'])
        max_l = str(summ['max_lambda'])
        mystr = (
            "initial_state=" + init_st + " " +
            "final_state=" + fin_st + " " +
            "min_lambda=" + min_l + " " +
            "max_lambda=" + max_l + " "
        )
        return mystr


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
        return AllInXEnsemble(volume, 0)

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
        return AllInXEnsemble(volume, -1)

    @staticmethod
    def A2BEnsemble(volume_a, volume_b, trusted = True):
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
            SingleFrameEnsemble(AllInXEnsemble(volume_a)),
            AllOutXEnsemble(volume_a | volume_b),
            SingleFrameEnsemble(AllInXEnsemble(volume_b))
        ])


    @staticmethod
    def TISEnsembleSet(volume_a, volume_b, volumes_x, orderparameter):
        myset = []
        for vol in volumes_x:
            myset.append(
                paths.TISEnsemble(volume_a, volume_b, vol, orderparameter)
            )
        return myset

