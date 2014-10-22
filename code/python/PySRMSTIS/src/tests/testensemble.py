from nose.tools import assert_equal, assert_not_equal, raises
from nose.plugins.skip import SkipTest
from duckpunching import CallIdentity, prepend_exception_message

import os
import sys
sys.path.append(os.path.abspath('../'))
from trajectory import Trajectory
from volume import LambdaVolume
from ensemble import *

import re
import random

def wrap_traj(traj, start, length):
    """Wraps the traj such that the original traj starts at frame `start`
    and is of length `length` by padding beginning with traj[0] and end with
    traj[-1]. Used to test the slice restricted trajectories."""
    if (start < 0) or (length < len(traj)+start):
        raise ValueError("""wrap_traj: start < 0 or length < len(traj)+start
                                {0} < 0 or {1} < {2}+{0}""".format(
                                start, length, len(traj)))
    outtraj = traj[:] # shallow copy
    # prepend
    for i in range(start):
        outtraj.insert(0, traj[0])
    # append
    for i in range(length - (len(traj)+start)):
        outtraj.append(traj[-1])
    return outtraj

def test_wrap_traj():
    """Testing wrap_traj (oh gods, the meta! a test for a test function!)"""
    intraj = [1, 2, 3]
    assert_equal(wrap_traj(intraj, 3, 6), [1, 1, 1, 1, 2, 3])
    assert_equal(wrap_traj(intraj, 3, 8), [1, 1, 1, 1, 2, 3, 3, 3])
    assert_equal(wrap_traj(intraj, 3, 8)[slice(3, 6)], intraj)

def build_trajdict(trajtypes, lower, upper):
    upperadddict = {'a' : 'in', 'b' : 'out', 'c' : 'cross', 'o' : 'hit'}
    loweradddict = {'a' : 'out', 'b' : 'in', 'c' : 'in', 'o' : 'hit'}
    lowersubdict = {'a' : 'in', 'b' : 'out', 'c' : 'cross', 'o' : 'hit'}
    uppersubdict = {'a' : 'out', 'b' : 'in', 'c' : 'in', 'o' : 'hit'}
    adjustdict = {'a' : (lambda x: -0.05*x), 'b' : (lambda x: 0.05*x),
                  'c' : (lambda x: 0.05*x + 0.15), 'o' : (lambda x: 0.0)}
    mydict = {}
    for mystr in trajtypes:
        upperaddkey = "upper"
        uppersubkey = "upper"
        loweraddkey = "lower"
        lowersubkey = "lower"
        delta = []
        for char in mystr:
            upperaddkey += "_"+upperadddict[char]
            loweraddkey += "_"+loweradddict[char]
            uppersubkey += "_"+uppersubdict[char]
            lowersubkey += "_"+lowersubdict[char]
            delta.append(adjustdict[char](random.randint(1, 4)))

        mydict[upperaddkey] = map(upper.__add__, delta)
        mydict[loweraddkey] = map(lower.__add__, delta)
        mydict[uppersubkey] = map(upper.__sub__, delta)
        mydict[lowersubkey] = map(lower.__sub__, delta)
    return mydict 

def setUp():
    ''' Setup for tests of classes in ensemble.py. '''
    #random.seed
    global lower, upper, op, vol1, vol2, ttraj
    lower = 0.1
    upper = 0.5
    op = CallIdentity()
    vol1 = LambdaVolume(op, lower, upper)
    vol2 = LambdaVolume(op, -0.1, 0.7)
    # we use the following codes to describe trajectories:
    # in : in the state
    # out : out of the state
    # hit : on the state border
    #
    # deltas of each letter from state edge:
    # a < 0 ; 0 < b < 0.2 ; c > 0.2; o = 0
    trajtypes = ["a", "o", "ab", "aob", "bob", "aba", "aaa", "abcba",
                 "abaa", "abba", "abaab", "ababa", "abbab",
                 "abaaba", "aobab", "abab", "abcbababcba"]
    ttraj = build_trajdict(trajtypes, lower, upper)
    # TODO: add support for empty traj to current tests (or perhaps add a
    # separate test for empty traj)
    #ttraj['empty'] = []

    # make the tests from lists into trajectories
    for test in ttraj.keys():
        ttraj[test] = Trajectory(ttraj[test])

def in_out_parser(testname):
    allowed_parts = ['in', 'out']
    parts = re.split("_", testname)
    res = []
    for part in parts:
        to_append = None
        if part in allowed_parts:
            to_append = part
        elif part == 'hit':
            if 'upper' in parts:
                to_append = 'in'
            elif 'lower' in parts:
                to_append = 'in'
        elif part == 'cross':
            to_append = 'out'
        if to_append != None:
            if res == []:
                res.append(to_append)
            elif to_append != res[-1]:
                res.append(to_append)
    return res

class EnsembleTest(object):
    def _single_test(self, ensemble_fcn, traj, res, failmsg):
        try:
            assert_equal(ensemble_fcn(traj), res)
        except AssertionError as e:
            prepend_exception_message(e, failmsg)
            raise


    def _run(self, results):
        """Actually run tests on the trajectory and the wrapped trajectory.

        Nearly all of the tests are just this simple. By creating custom error
        messages (using prepend_exception_message) we can wrap the many tests
        into loops instead of making tons of lines of code.
        """
        for test in results.keys():
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.ensemble, ttraj[test], results[test], failmsg)

            wrapped = wrap_traj(ttraj[test], self.wrapstart, self.wrapend)
            lentt = len(ttraj[test])

            failmsg = "Failure in wrapped "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.ensemble, wrapped, results[test], failmsg)

            failmsg = "Failure in slice_ens "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.slice_ens, wrapped, results[test], failmsg)


class testLeaveXEnsemble(EnsembleTest):
    def setUp(self):
        self.leaveX = LeaveXEnsemble(vol1)

    def test_leaveX(self):
        """LeaveXEnsemble passes the trajectory test suite"""
        for test in ttraj.keys():
            if "out" in in_out_parser(test):
                res = True
            else:
                res = False
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.leaveX, ttraj[test], res, failmsg)

    def test_leaveX_str(self):
        volstr = "{x|Id(x) in [0.1, 0.5]}"
        assert_equal(self.leaveX.__str__(), 
                     "x[t] in (not "+volstr+") for one t in [1:-1]")

class testInXEnsemble(EnsembleTest):
    def setUp(self):
        self.inX = InXEnsemble(vol1)

    def test_inX(self):
        """InXEnsemble passes the trajectory test suite"""
        for test in ttraj.keys():
            if "out" in in_out_parser(test):
                res = False
            else:
                res = True
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.inX, ttraj[test], res, failmsg)

    def test_inX_str(self):
        volstr = "{x|Id(x) in [0.1, 0.5]}"
        assert_equal(self.inX.__str__(),
                     "x[t] in "+volstr+" for all t in [1:-1]")

class testOutXEnsemble(EnsembleTest):
    def setUp(self):
        self.outX = OutXEnsemble(vol1)

    def test_outX(self):
        """OutXEnsemble passes the trajectory test suite"""
        for test in ttraj.keys():
            if "in" in in_out_parser(test):
                res = False
            else:
                res = True
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.outX, ttraj[test], res, failmsg)

    def test_outX_str(self):
        volstr = "{x|Id(x) in [0.1, 0.5]}"
        assert_equal(self.outX.__str__(),
                     "x[t] in (not "+volstr+") for all t in [1:-1]")

class testHitXEnsemble(EnsembleTest):
    def setUp(self):
        self.hitX = HitXEnsemble(vol1)

    def test_hitX(self):
        """HitXEnsemble passes the trajectory test suite"""
        for test in ttraj.keys():
            if "in" in in_out_parser(test):
                res = True
            else:
                res = False
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.hitX, ttraj[test], res, failmsg)

    def test_hitX_str(self):
        volstr = "{x|Id(x) in [0.1, 0.5]}"
        assert_equal(self.hitX.__str__(),
                     "x[t] in "+volstr+" for one t in [1:-1]")

class testExitsXEnsemble(EnsembleTest):
    def setUp(self):
        self.ensemble = ExitsXEnsemble(vol1)
        # longest ttraj is 6 = 9-3 frames long
        self.slice_ens = ExitsXEnsemble(vol1, slice(3,9))
        self.wrapstart = 3
        self.wrapend = 12

    @raises(ValueError)
    def test_single_frame(self):
        '''ExitsXEnsemble built with single frame raises ValueError'''
        single_frame_ensemble = ExitsXEnsemble(vol1, 2)

    @raises(ValueError)
    def test_single_frame_by_modification(self):
        '''ExitsXEnsemble modified to have single frame raises ValueError'''
        single_frame_ensemble = ExitsXEnsemble(vol1)
        single_frame_ensemble.frames = 2
        single_frame_ensemble(ttraj['upper_in_out'])

    def test_noncrossing(self):
        '''ExitsXEnsemble for noncrossing trajectories'''
        results = { 'upper_in' : False,
                    'upper_out' : False,
                    'lower_in' : False,
                    'lower_out' : False
                  }
        self._run(results)

    def test_hitsborder(self):
        '''ExitsXEnsemble for border-hitting trajectories'''
        results = { 'lower_in_hit_in' : False,
                    'upper_in_hit_in' : False,
                    'lower_out_hit_out' : True,
                    'upper_out_hit_out' : True
                  }
        self._run(results)

    def test_exit(self):
        '''ExitsXEnsemble for exiting trajecories'''
        results = { 'lower_in_out' : True,
                    'upper_in_out' : True,
                    'lower_in_hit_out' : True,
                    'upper_in_hit_out' : True,
                    'lower_out_in_out_in' : True,
                    'upper_out_in_out_in' : True,
                    'lower_in_out_in_out' : True,
                    'upper_in_out_in_out' : True
                  }
        self._run(results)

    def test_entrance(self):
        '''ExitsXEnsemble for entering trajectories'''
        results = { 'lower_out_in' : False,
                    'upper_out_in' : False,
                    'lower_out_hit_in' : False,
                    'upper_out_hit_in' : False,
                    'lower_out_in_out_in' : True,
                    'upper_out_in_out_in' : True,
                    'lower_in_out_in_out' : True,
                    'upper_in_out_in_out' : True
                  }
        self._run(results)

    def test_str(self):
        assert_equal(self.ensemble.__str__(),
            'exists x[t], x[t+1] in [None:None] such that x[t] in {0} and x[t+1] not in {0}'.format(vol1))
        self.ensemble.frames = slice(3,8)
        assert_equal(self.ensemble.__str__(),
            'exists x[t], x[t+1] in [3:8] such that x[t] in {0} and x[t+1] not in {0}'.format(vol1))

class testEntersXEnsemble(testExitsXEnsemble):
    def setUp(self):
        self.ensemble = EntersXEnsemble(vol1)
        # longest ttraj is 6 = 9-3 frames long
        self.slice_ens = EntersXEnsemble(vol1, slice(3,9))
        self.wrapstart = 3
        self.wrapend = 12

    @raises(ValueError)
    def test_single_frame(self):
        '''EntersXEnsemble built with single frame raises ValueError'''
        single_frame_ensemble = EntersXEnsemble(vol1, 2)

    @raises(ValueError)
    def test_single_frame_by_modification(self):
        '''EntersXEnsemble modified to have single frame raises ValueError'''
        single_frame_ensemble = EntersXEnsemble(vol1)
        single_frame_ensemble.frames = 2
        single_frame_ensemble(ttraj['upper_in_out'])

    def test_noncrossing(self):
        '''EntersXEnsemble for noncrossing trajectories'''
        results = { 'upper_in_in_in' : False,
                    'upper_out_out_out' : False,
                    'lower_in_in_in' : False,
                    'lower_out_out_out' : False
                  }
        self._run(results)

    def test_hitsborder(self):
        '''EntersXEnsemble for border-hitting trajectories'''
        results = { 'lower_in_hit_in' : False,
                    'upper_in_hit_in' : False,
                    'lower_out_hit_out' : True,
                    'upper_out_hit_out' : True
                  }
        self._run(results)

    def test_exit(self):
        '''EntersXEnsemble for exiting trajecories'''
        results = { 'lower_in_out' : False,
                    'upper_in_out' : False,
                    'lower_in_hit_out' : False,
                    'upper_in_hit_out' : False,
                    'lower_out_in_out_in' : True,
                    'upper_out_in_out_in' : True,
                    'lower_in_out_in_out' : True,
                    'upper_in_out_in_out' : True
                  }
        self._run(results)

    def test_entrance(self):
        '''EntersXEnsemble for entering trajectories'''
        results = { 'lower_out_in' : True,
                    'upper_out_in' : True,
                    'lower_out_hit_in' : True,
                    'upper_out_hit_in' : True,
                    'lower_out_in_out_in' : True,
                    'upper_out_in_out_in' : True,
                    'lower_in_out_in_out' : True,
                    'upper_in_out_in_out' : True
                  }
        self._run(results)

    def test_str(self):
        assert_equal(self.ensemble.__str__(),
            'exists x[t], x[t+1] in [None:None] such that x[t] not in {0} and x[t+1] in {0}'.format(vol1))
        self.ensemble.frames = slice(3,8)
        assert_equal(self.ensemble.__str__(),
            'exists x[t], x[t+1] in [3:8] such that x[t] not in {0} and x[t+1] in {0}'.format(vol1))

class testSequentialEnsembles(EnsembleTest):
    def setUp(self):
        self.inX = InXEnsemble(vol1, frames=slice(None,None))
        self.outX = OutXEnsemble(vol1, frames=slice(None,None))
        self.hitX = HitXEnsemble(vol1, frames=slice(None,None))
        self.leaveX = LeaveXEnsemble(vol1, frames=slice(None,None))
        self.enterX = EntersXEnsemble(vol1, frames=slice(None,None))
        self.exitX = ExitsXEnsemble(vol1, frames=slice(None,None))
        self.inInterface = InXEnsemble(vol2, frames=slice(None,None))
        self.leaveX0 = LeaveXEnsemble(vol2, frames=slice(None,None))
        self.inX0 = InXEnsemble(vol2, frames=slice(None,None))
        self.length1 = LengthEnsemble(1)
        self.length0 = LengthEnsemble(0)
        # pseudo_tis and pseudo_minus assume that the interface is equal to
        # the state boundary
        self.pseudo_tis = SequentialEnsemble( [
                                    self.inX & self.length1,
                                    self.outX,
                                    self.inX & self.length1 ]
                                    )
        self.pseudo_minus = SequentialEnsemble( [
                                    self.inX & self.length1,
                                    self.outX,
                                    self.inX,
                                    self.outX,
                                    self.inX & self.length1 ]
        )
        self.minus = SequentialEnsemble([
            self.inX & self.length1,
            self.outX & self.leaveX0,
            self.inX & self.length1,
            self.inX0 | self.length0,
            self.outX & self.leaveX0,
            self.inX & self.length1
                                    ])
        self.tis = SequentialEnsemble([
            self.inX & self.length1,
            self.outX & self.leaveX0,
            self.inX & self.length1,
        ])

    @raises(ValueError)
    def test_maxminoverlap_size(self):
        """SequentialEnsemble errors if max/min overlap sizes different"""
        SequentialEnsemble([self.inX, self.outX, self.inX], (0,0), (0,0,0))

    @raises(ValueError)
    def test_maxoverlap_ensemble_size(self):
        """SequentialEnsemble errors if overlap sizes don't match ensemble size"""
        SequentialEnsemble([self.inX, self.outX, self.inX], (0,0,0), (0,0,0))

    @raises(ValueError)
    def test_minmax_order(self):
        """SequentialEnsemble errors if min_overlap > max_overlap"""
        SequentialEnsemble([self.inX, self.outX, self.inX], (0,1), (0,0))

    def test_allowed_initializations(self):
        """SequentialEnsemble initializes correctly with defaults"""
        A = SequentialEnsemble([self.inX, self.outX, self.inX], (0,0), (0,0))
        B = SequentialEnsemble([self.inX, self.outX, self.inX],0,0)
        C = SequentialEnsemble([self.inX, self.outX, self.inX])
        assert_equal(A.min_overlap,B.min_overlap)
        assert_equal(A.min_overlap,C.min_overlap)
        assert_equal(A.max_overlap,B.max_overlap)
        assert_equal(A.max_overlap,C.max_overlap)

    def test_overlap_max(self):
        """SequentialEnsemble allows overlaps up to overlap max, no more"""
        raise SkipTest

    def test_overlap_min(self):
        """SequentialEnsemble requires overlaps of at least overlap min"""
        raise SkipTest

    def test_overlap_max_inf(self):
        """SequentialEnsemble works if max overlap in infinite"""
        raise SkipTest

    def test_overlap_min_gap(self):
        """SequentialEnsemble works in mix overlap is negative (gap)"""
        raise SkipTest

    def test_overlap_max_gap(self):
        """SequentialEnsemble works if max overlap is negative (gap)"""
        raise SkipTest

    def test_can_append_tis(self):
        """SequentialEnsemble as TISEnsemble knows when it can append"""
        results =   {   'upper_in_out' : True,
                        'lower_in_out' : True,
                        'upper_in_out_in' : False,
                        'lower_in_out_in' : False,
                        'upper_in' : True,
                        'lower_in' : True,
                        'upper_in_in_in' : False,
                        'lower_in_in_in' : False,
                        'upper_out_out_out' : True,
                        'lower_out_out_out' : True,
                        'upper_out_in' : False,
                        'lower_out_in' : False,
                        'upper_out' : True,
                        'lower_out' : True,
                        'upper_in_out_in_in' : False,
                        'lower_in_out_in_in' : False,
                        'upper_in_out_in_out_in' : False,
                        'lower_in_out_in_out_in' : False
                    }   
        for test in results.keys():
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.pseudo_tis.forward, 
                                ttraj[test], results[test], failmsg)

    def test_can_append_minus(self):
        """SequentialEnsemble as MinusEnsemble knows when it can append"""
        results =   {   'upper_in_out' : True,
                        'lower_in_out' : True,
                        'upper_in_out_in' : True,
                        'lower_in_out_in' : True,
                        'upper_in' : True,
                        'lower_in' : True,
                        'upper_in_in_in' : True,
                        'lower_in_in_in' : True,
                        'upper_out_out_out' : True,
                        'lower_out_out_out' : True,
                        'upper_out_in' : True,
                        'lower_out_in' : True,
                        'upper_out' : True,
                        'lower_out' : True,

                        'upper_in_out_in_in' : True,
                        'lower_in_out_in_in' : True,
                        'upper_in_out_in_out_in' : False,
                        'lower_in_out_in_out_in' : False,
                        'upper_in_out_in_in_out' : True,
                        'lower_in_out_in_in_out' : True,
                        'upper_out_in_out' : True,
                        'lower_out_in_out' : True,
                        'upper_out_in_in_out' : True,
                        'lower_out_in_in_out' : True,
                        'upper_out_in_out_in': False,
                        'lower_out_in_out_in': False,
                        'upper_out_in_in_out_in' : False,
                        'lower_out_in_in_out_in' : False
                    }   
        for test in results.keys():
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.pseudo_minus.forward, 
                                ttraj[test], results[test], failmsg)

    def test_can_prepend_tis(self):
        """SequentialEnsemble as TISEnsemble knows when it can prepend"""
        results =   {   'upper_in_out' : False,
                        'lower_in_out' : False,
                        'upper_in_out_in' : False,
                        'lower_in_out_in' : False,
                        'upper_in' : True,
                        'lower_in' : True,
                        'upper_in_in_in' : False,
                        'lower_in_in_in' : False,
                        'upper_out_out_out' : True,
                        'lower_out_out_out' : True,
                        'upper_out_in' : True,
                        'lower_out_in' : True,
                        'upper_out' : True,
                        'lower_out' : True,
                        'upper_in_out_in_in' : False,
                        'lower_in_out_in_in' : False,
                        'upper_in_out_in_out_in' : False,
                        'lower_in_out_in_out_in' : False
                    }   
        for test in results.keys():
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.pseudo_tis.backward, 
                                ttraj[test], results[test], failmsg)


    def test_can_prepend_minus(self):
        """SequentialEnsemble as MinusEnsemble knows when it can prepend"""
        results =   {   'upper_in_out' : True,
                        'lower_in_out' : True,
                        'upper_in_out_in' : True,
                        'lower_in_out_in' : True,
                        'upper_in' : True,
                        'lower_in' : True,
                        'upper_in_in_in' : True,
                        'lower_in_in_in' : True,
                        'upper_out_out_out' : True,
                        'lower_out_out_out' : True,
                        'upper_out_in' : True,
                        'lower_out_in' : True,
                        'upper_out' : True,
                        'lower_out' : True,

                        'upper_in_out_in_in' : False,
                        'lower_in_out_in_in' : False,
                        'upper_in_out_in_out_in' : False,
                        'lower_in_out_in_out_in' : False,
                        'upper_in_out_in_in_out' : False,
                        'lower_in_out_in_in_out' : False,
                        'upper_out_in_out' : True,
                        'lower_out_in_out' : True,
                        'upper_out_in_in_out' : True,
                        'lower_out_in_in_out' : True,
                        'upper_out_in_out_in': True,
                        'lower_out_in_out_in': True,
                        'upper_out_in_in_out_in' : True,
                        'lower_out_in_in_out_in' : True
                    }   
        for test in results.keys():
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.pseudo_minus.backward, 
                                ttraj[test], results[test], failmsg)

    
    def test_sequential_transition_frames(self):
        """SequentialEnsemble identifies transitions frames correctly"""
        ensemble = SequentialEnsemble([self.inX, self.outX])
        results = {'upper_in_in_in' : [3],
                   'upper_out_out_out' : [],
                   'upper_in_out_in' : [1,2],
                   'upper_in_out' : [1,2]
                  }
        for test in results.keys():
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(ensemble.transition_frames, 
                                ttraj[test], results[test], failmsg)

    def test_sequential_simple_in_out_call(self):
        """Simplest sequential ensemble identifies correctly"""
        ensemble = SequentialEnsemble([self.inX, self.outX])
        results = {'upper_in_in_in' : False,
                   'upper_out_out_out' : False,
                   'upper_in_out_in' : False,
                   'upper_in_out' : True
                  }
        print
        for test in results.keys():
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            print test,
            self._single_test(ensemble, 
                                ttraj[test], results[test], failmsg)



    def test_sequential_in_out(self):
        """SequentialEnsembles based on In/OutXEnsemble"""
        # idea: for each ttraj, use the key name to define in/out behavior,
        # dynamically construct a SequentialEnsemble
        ens_dict = {'in' : self.inX, 'out' : self.outX }
        for test in ttraj.keys():
            ens_list = in_out_parser(test)
            ens = []

            # how to pick ensembles is specific to this test
            for ens_type in ens_list:
                ens.append(ens_dict[ens_type])

            ensemble = SequentialEnsemble(ens)
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(ensemble, ttraj[test], True, failmsg)

    def test_sequential_pseudo_tis(self):
        """SequentialEnsemble as Pseudo-TISEnsemble identifies paths"""
        results = {}
        for test in ttraj.keys():
            results[test] = False
        results['upper_in_out_in'] = True
        results['lower_in_out_in'] = True
        results['upper_in_out_out_in'] = True
        results['lower_in_out_out_in'] = True
        results['lower_in_out_cross_out_in'] = True
        results['upper_in_out_cross_out_in'] = True
        for test in results.keys():
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.pseudo_tis, ttraj[test], results[test],
                              failmsg)

    def test_sequential_pseudo_minus(self):
        """SequentialEnsemble as Pseudo-MinusEnsemble identifies paths"""
        results = {}
        for test in ttraj.keys():
            results[test] = False
        results['upper_in_out_in_out_in'] = True
        results['lower_in_out_in_out_in'] = True
        results['upper_in_out_in_in_out_in'] = True
        results['lower_in_out_in_in_out_in'] = True
        for test in results.keys():
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.pseudo_minus, ttraj[test], results[test],
                              failmsg)

    def test_sequential_tis(self):
        """SequentialEnsemble as TISEnsemble identifies paths"""
        raise SkipTest
        results = {}
        for test in ttraj.keys():
            results[test] = False
        results['upper_in_out_cross_out_in'] = True
        results['lower_in_out_cross_out_in'] = True
        for test in results.keys():
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.tis, ttraj[test], results[test],
                              failmsg)


    def test_sequential_enter_exit(self):
        """SequentialEnsembles based on Enters/ExitsXEnsemble"""
        # TODO: this includes a test of the overlap ability
        raise SkipTest


    def test_str(self):
        assert_equal(self.pseudo_tis.__str__(), """[
(
  x[t] in {x|Id(x) in [0.1, 0.5]} for all t in [None:None]
)
and
(
  len(x) = 1
),
x[t] in (not {x|Id(x) in [0.1, 0.5]}) for all t in [None:None],
(
  x[t] in {x|Id(x) in [0.1, 0.5]} for all t in [None:None]
)
and
(
  len(x) = 1
)
]""")
        self.pseudo_minus.__str__()
