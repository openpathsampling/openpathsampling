from nose.tools import assert_equal, assert_not_equal, raises
from nose.plugins.skip import SkipTest
from test_helpers import CallIdentity, prepend_exception_message

import openpathsampling as paths
from openpathsampling.ensemble import *

import logging
logging.getLogger('opentis.ensemble').setLevel(logging.DEBUG)
logging.getLogger('opentis.initialization').setLevel(logging.CRITICAL)

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
                  'c' : (lambda x: 0.05*x + 0.16), 'o' : (lambda x: 0.0)}
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

def tstr(ttraj):
    return list(ttraj).__str__()

def results_upper_lower(adict):
    res_dict = {}
    for test in adict.keys():
        res_dict['upper_'+test] = adict[test]
        res_dict['lower_'+test] = adict[test]
    return res_dict


def setUp():
    ''' Setup for tests of classes in ensemble.py. '''
    #random.seed
    global lower, upper, op, vol1, vol2, vol3, ttraj
    lower = 0.1
    upper = 0.5
    op = CallIdentity()
    vol1 = paths.LambdaVolume(op, lower, upper)
    vol2 = paths.LambdaVolume(op, -0.1, 0.7)
    vol3 = paths.LambdaVolume(op, 2.0, 2.5)
    # we use the following codes to describe trajectories:
    # in : in the state
    # out : out of the state
    # hit : on the state border
    #
    # deltas of each letter from state edge:
    # a < 0 ; 0 < b < 0.2 ; c > 0.2; o = 0
    trajtypes = ["a", "o", "aa", "ab", "aob", "bob", "aba", "aaa", "abcba",
                 "abaa", "abba", "abaab", "ababa", "abbab", "ac", "bc",
                 "abaaba", "aobab", "abab", "abcbababcba", "aca", "abc",
                 "acaca", "acac", "caca", "aaca", "baca", "aaba", "aab",
                 "aabbaa"
                ]
    ttraj = build_trajdict(trajtypes, lower, upper)

    # make the tests from lists into trajectories
    for test in ttraj.keys():
        ttraj[test] = paths.Trajectory(ttraj[test])

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
    def __init__(self):
        self.length0 = LengthEnsemble(0)

    def _single_test(self, ensemble_fcn, traj, res, failmsg):
        try:
            assert_equal(ensemble_fcn(traj), res)
        except AssertionError as e:
            prepend_exception_message(e, failmsg)
            raise

    def _test_everything(self, test_fcn, non_default=[], default=False):
        """
        Runs tests using *all* the trajectory test suite. This is the
        ultimate in test-running simplicity!!
        """
        results = {}
        for test in ttraj.keys():
            results[test] = default
        nondef_dict = {}
        for test in non_default:
            if test in ttraj.keys():
                results[test] = not default
            if "lower_"+test in ttraj.keys():
                results["lower_"+test] = not default
            if "upper_"+test in ttraj.keys():
                results["upper_"+test] = not default

        for test in results.keys():
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(test_fcn, ttraj[test], results[test], failmsg)


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

    def test_leaveX_0(self):
        """LeaveXEnsemble treatment of zero-length trajectory"""
        assert_equal(self.leaveX(paths.Trajectory([])), False)
        assert_equal(self.leaveX.can_append(paths.Trajectory([])), True)
        assert_equal(self.leaveX.can_prepend(paths.Trajectory([])), True)

    def test_leaveX_str(self):
        volstr = "{x|Id(x) in [0.1, 0.5]}"
        assert_equal(self.leaveX.__str__(), 
                     "exists t such that x[t] in (not "+volstr+")")

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

    def test_inX_0(self):
        """InXEnsemble treatment of zero-length trajectory"""
        assert_equal(self.inX(paths.Trajectory([])), False)
        assert_equal(self.inX.can_append(paths.Trajectory([])), True)
        assert_equal(self.inX.can_prepend(paths.Trajectory([])), True)

    def test_inX_str(self):
        volstr = "{x|Id(x) in [0.1, 0.5]}"
        assert_equal(self.inX.__str__(),
                     "x[t] in "+volstr+" for all t")

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

    def test_outX_0(self):
        """OutXEnsemble treatment of zero-length trajectory"""
        assert_equal(self.outX(paths.Trajectory([])), False)
        assert_equal(self.outX.can_append(paths.Trajectory([])), True)
        assert_equal(self.outX.can_prepend(paths.Trajectory([])), True)

    def test_outX_str(self):
        volstr = "{x|Id(x) in [0.1, 0.5]}"
        assert_equal(self.outX.__str__(),
                     "x[t] in (not "+volstr+") for all t")

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

    def test_hitX_0(self):
        """HitXEnsemble treatment of zero-length trajectory"""
        assert_equal(self.hitX(paths.Trajectory([])), False)
        assert_equal(self.hitX.can_append(paths.Trajectory([])), True)
        assert_equal(self.hitX.can_prepend(paths.Trajectory([])), True)

    def test_hitX_str(self):
        volstr = "{x|Id(x) in [0.1, 0.5]}"
        assert_equal(self.hitX.__str__(),
                     "exists t such that x[t] in "+volstr)

class testExitsXEnsemble(EnsembleTest):
    def setUp(self):
        self.ensemble = ExitsXEnsemble(vol1)
        # longest ttraj is 6 = 9-3 frames long
        self.slice_ens = ExitsXEnsemble(vol1, slice(3,9))
        self.wrapstart = 3
        self.wrapend = 12

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
            'exists x[t], x[t+1] such that x[t] in {0} and x[t+1] not in {0}'.format(vol1))

class testEntersXEnsemble(testExitsXEnsemble):
    def setUp(self):
        self.ensemble = EntersXEnsemble(vol1)
        # longest ttraj is 6 = 9-3 frames long
        self.slice_ens = EntersXEnsemble(vol1, slice(3,9))
        self.wrapstart = 3
        self.wrapend = 12

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
            'exists x[t], x[t+1] such that x[t] not in {0} and x[t+1] in {0}'.format(vol1))

class testSequentialEnsemble(EnsembleTest):
    def setUp(self):
        self.inX = InXEnsemble(vol1)
        self.outX = OutXEnsemble(vol1)
        self.hitX = HitXEnsemble(vol1)
        self.leaveX = LeaveXEnsemble(vol1)
        self.enterX = EntersXEnsemble(vol1)
        self.exitX = ExitsXEnsemble(vol1)
        self.inInterface = InXEnsemble(vol2)
        self.leaveX0 = LeaveXEnsemble(vol2)
        self.inX0 = InXEnsemble(vol2)
        self.length1 = LengthEnsemble(1)
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
            self._single_test(self.pseudo_tis.can_append, 
                                ttraj[test], results[test], failmsg)

    def test_can_append_pseudominus(self):
        """SequentialEnsemble as Pseudo-MinusEnsemble knows when it can append"""
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
                        'lower_out_in_in_out_in' : False,
                        'upper_in_cross_in' : True,
                        'lower_in_cross_in' : True,
                        'upper_in_cross_in_cross' : True,
                        'lower_in_cross_in_cross' : True,
                        'upper_cross_in_cross_in' : False,
                        'lower_cross_in_cross_in' : False,
                        'upper_in_cross_in_cross_in' : False,
                        'lower_in_cross_in_cross_in' : False
                    }   
        for test in results.keys():
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.pseudo_minus.can_append, 
                                ttraj[test], results[test], failmsg)

    def test_can_append_minus(self):
        results = {'upper_in_cross_in' : True,
                   'lower_in_cross_in' : True,
                   'upper_in_cross_in_cross' : True,
                   'lower_in_cross_in_cross' : True,
                   'upper_cross_in_cross_in' : False,
                   'lower_cross_in_cross_in' : False,
                   'upper_in_cross_in_cross_in' : False,
                   'lower_in_cross_in_cross_in' : False
                  }
        for test in results.keys():
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.minus.can_append, 
                                ttraj[test], results[test], failmsg)

    def test_can_prepend_minus(self):
        results = {'upper_in_cross_in' : True,
                   'lower_in_cross_in' : True,
                   'upper_in_cross_in_cross' : False,
                   'lower_in_cross_in_cross' : False,
                   'upper_cross_in_cross_in' : True,
                   'lower_cross_in_cross_in' : True,
                   'upper_in_cross_in_cross_in' : False,
                   'lower_in_cross_in_cross_in' : False
                  }
        for test in results.keys():
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.minus.can_append, 
                                ttraj[test], results[test], failmsg)

    def test_can_prepend_pseudo_tis(self):
        """SequentialEnsemble as Pseudo-TISEnsemble knows when it can prepend"""
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
            self._single_test(self.pseudo_tis.can_prepend, 
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
            self._single_test(self.pseudo_minus.can_prepend, 
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
        for test in results.keys():
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
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
        results['upper_in_cross_in'] = True
        results['lower_in_cross_in'] = True
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
        results['upper_in_cross_in_cross_in'] = True
        results['lower_in_cross_in_cross_in'] = True
        for test in results.keys():
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.pseudo_minus, ttraj[test], results[test],
                              failmsg)

    def test_sequential_tis(self):
        """SequentialEnsemble as TISEnsemble identifies paths"""
        results = {}
        for test in ttraj.keys():
            results[test] = False
        results['upper_in_out_cross_out_in'] = True
        results['lower_in_out_cross_out_in'] = True
        results['upper_in_cross_in'] = True
        results['lower_in_cross_in'] = True
        for test in results.keys():
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.tis, ttraj[test], results[test], failmsg)
    
    def test_sequential_minus(self):
        """SequentialEnsemble as MinusEnsemble identifies paths"""
        results = {}
        for test in ttraj.keys():
            results[test] = False
        results['upper_in_out_cross_out_in_out_in_out_cross_out_in'] = True
        results['lower_in_out_cross_out_in_out_in_out_cross_out_in'] = True
        results['upper_in_cross_in_cross_in'] = True
        results['lower_in_cross_in_cross_in'] = True

        for test in results.keys():
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.minus, ttraj[test], results[test], failmsg)

    def test_sequential_generate_first_tis(self):
        """SequentialEnsemble to generate the first TIS path"""
        ensemble = SequentialEnsemble([
            self.outX | self.length0,
            self.inX,
            self.outX & self.leaveX0,
            self.inX & self.length1
        ])
        match_results = {
            'upper_in_in_cross_in' : True,
            'lower_in_in_cross_in' : True,
            'upper_out_in_cross_in' : True,
            'lower_out_in_cross_in' : True,
            'upper_in_cross_in' : True,
            'lower_in_cross_in' : True
        }
        logging.getLogger('opentis.ensemble').info("Starting tests....")
        for test in match_results.keys():
            failmsg = "Match failure in "+test+"("+str(ttraj[test])+"): "
            logging.getLogger('opentis.ensemble').info(
                "Testing: "+str(test)
            )
            self._single_test(ensemble, ttraj[test], 
                              match_results[test], failmsg)

        append_results = {
            'upper_in' : True,
            'upper_in_in_in' : True,
            'upper_in_in_out' : True,
            'upper_in_in_out_in' : False,
            'upper_out' : True,
            'upper_out_out_out' : True,
            'upper_out_in_out' : True,
            'upper_out_in_out_in' : False
        }
        for test in append_results.keys():
            failmsg = "Append failure in "+test+"("+str(ttraj[test])+"): "
            logging.getLogger('opentis.ensemble').info(
                "Testing: "+str(test)
            )
            self._single_test(ensemble.can_append, ttraj[test], 
                              append_results[test], failmsg)

    def test_sequential_enter_exit(self):
        """SequentialEnsembles based on Enters/ExitsXEnsemble"""
        # TODO: this includes a test of the overlap ability
        raise SkipTest


    def test_str(self):
        assert_equal(self.pseudo_tis.__str__(), """[
(
  x[t] in {x|Id(x) in [0.1, 0.5]} for all t
)
and
(
  len(x) = 1
),
x[t] in (not {x|Id(x) in [0.1, 0.5]}) for all t,
(
  x[t] in {x|Id(x) in [0.1, 0.5]} for all t
)
and
(
  len(x) = 1
)
]""")

class testSlicedTrajectoryEnsemble(EnsembleTest):
    def test_sliced_ensemble_init(self):
        init_as_int = SlicedTrajectoryEnsemble(InXEnsemble(vol1), 3)
        init_as_slice = SlicedTrajectoryEnsemble(InXEnsemble(vol1),
                                                 slice(3, 4))
        assert_equal(init_as_int, init_as_slice)
        assert_equal(init_as_slice.slice, init_as_int.slice)

    def test_sliced_as_TISEnsemble(self):
        '''SlicedTrajectory and Sequential give same TIS results'''
        sliced_tis = (
            SlicedTrajectoryEnsemble(InXEnsemble(vol1), 0) &
            SlicedTrajectoryEnsemble(OutXEnsemble(vol1 | vol3), slice(1,-1)) & 
            SlicedTrajectoryEnsemble(LeaveXEnsemble(vol2), slice(1,-1)) &
            SlicedTrajectoryEnsemble(InXEnsemble(vol1 | vol3), -1)
        )
        sequential_tis = SequentialEnsemble([
            InXEnsemble(vol1) & LengthEnsemble(1),
            OutXEnsemble(vol1 | vol3) & LeaveXEnsemble(vol2),
            InXEnsemble(vol1 | vol3) & LengthEnsemble(1)
        ])
        for test in ttraj.keys():
            failmsg = "Failure in "+test+"("+tstr(ttraj[test])+"): "
            self._single_test(sliced_tis, ttraj[test], 
                              sequential_tis(ttraj[test]), failmsg)

    def test_slice_outside_trajectory_range(self):
        ens = SlicedTrajectoryEnsemble(InXEnsemble(vol1), slice(5,9))
        test = 'upper_in'
        # the slice should return the empty trajectory, and therefore should
        # return false
        assert_equal(ens(ttraj[test]), False)

    def test_even_sliced_trajectory(self):
        even_slice = slice(None, None, 2)
        ens = SlicedTrajectoryEnsemble(InXEnsemble(vol1), even_slice)
        bare_results = {'in' : True,
                        'in_in' : True,
                        'in_in_in' : True,
                        'in_out_in' : True,
                        'in_in_out' : False,
                        'in_out_in_in' : True,
                        'in_out_in_out_in' : True,
                        'out' : False,
                        'in_in_out_in' : False,
                        'in_cross_in_cross' : True
                       }
        results = results_upper_lower(bare_results)
        for test in results.keys():
            failmsg = "Failure in "+test+"("+tstr(ttraj[test])+"): "
            self._single_test(ens, ttraj[test], results[test], failmsg)

    def test_sliced_sequential_global_whole(self):
        even_slice = slice(None, None, 2)
        ens = SlicedTrajectoryEnsemble(SequentialEnsemble([
            InXEnsemble(vol1),
            OutXEnsemble(vol1)
        ]), even_slice)

        bare_results = {'in_in_out' : True,
                        'in_hit_out' : True,
                        'in_out_out_in_out' : True,
                        'in_hit_out_in_out' : True,
                        'in_out_out_in' : True,
                        'in_cross_in_cross' : False,
                        'in_out_out_in' : True,
                        'in_out_in' : False
                       }
        results = results_upper_lower(bare_results)
        for test in results.keys():
            failmsg = "Failure in "+test+"("+tstr(ttraj[test])+"): "
            self._single_test(ens, ttraj[test], results[test], failmsg)

    def test_sliced_sequential_subtraj_member(self):
        even_slice = slice(None, None, 2)
        ens = SequentialEnsemble([
            InXEnsemble(vol1),
            SlicedTrajectoryEnsemble(OutXEnsemble(vol1), even_slice)
        ])
        bare_results = {'in_out_in' : True,
                        'in_out_out_in' : False,
                        'in_in_out_in' : True,
                        'in_in_out' : True,
                        'in_in_cross_in' : True,
                        'in_out_in_out' : True,
                        'in_out_cross_out_in_out_in_out_cross_out_in' : True,
                        'in_out_in_in_out_in' : False
                       }
        results = results_upper_lower(bare_results)
        for test in results.keys():
            failmsg = "Failure in "+test+"("+tstr(ttraj[test])+"): "
            self._single_test(ens, ttraj[test], results[test], failmsg)

    def test_sliced_sequential_subtraj_middle(self):
        even_slice = slice(None, None, 2)
        ens = SequentialEnsemble([
            InXEnsemble(vol1),
            SlicedTrajectoryEnsemble(OutXEnsemble(vol1), even_slice),
            InXEnsemble(vol1) & LengthEnsemble(1)
        ])
        bare_results = {'in_in_out_out_in_in' : False
                       }
        results = results_upper_lower(bare_results)
        for test in results.keys():
            failmsg = "Failure in "+test+"("+tstr(ttraj[test])+"): "
            self._single_test(ens, ttraj[test], results[test], failmsg)


    def test_sliced_str(self):
        even_slice = slice(None,None, 2)
        slice_1_10 = slice(1, 10)
        slice_1_end = slice(1,None)
        slice_no_ends = slice(1, -1)
        inX = InXEnsemble(vol1)
        inXstr = "x[t] in {x|Id(x) in [0.1, 0.5]} for all t"
        assert_equal(SlicedTrajectoryEnsemble(inX, even_slice).__str__(),
                     "("+inXstr+" in {:} every 2)")
        assert_equal(SlicedTrajectoryEnsemble(inX, slice_1_10).__str__(),
                     "("+inXstr+" in {1:10})")
        assert_equal(SlicedTrajectoryEnsemble(inX, slice_1_end).__str__(),
                     "("+inXstr+" in {1:})")
        assert_equal(SlicedTrajectoryEnsemble(inX, slice_no_ends).__str__(),
                     "("+inXstr+" in {1:-1})")

class testOptionalEnsemble(EnsembleTest):
    def setUp(self):
        self.start_opt = SequentialEnsemble([
            OptionalEnsemble(OutXEnsemble(vol1)),
            InXEnsemble(vol1),
            OutXEnsemble(vol1),
        ])
        self.end_opt = SequentialEnsemble([
            OutXEnsemble(vol1),
            InXEnsemble(vol1),
            OptionalEnsemble(OutXEnsemble(vol1))
        ])
        self.mid_opt = SequentialEnsemble([
            InXEnsemble(vol1),
            OptionalEnsemble(OutXEnsemble(vol1) & InXEnsemble(vol2)),
            OutXEnsemble(vol2),
        ])

    def test_optional_start(self):
        bare_results = {'in_out' : True,
                        'in_in_out' : True,
                        'out_in_out' : True,
                        'out_out' : False,
                        'out_in_in_out' : True,
                        'in_out_in' : False
                       }
        results = results_upper_lower(bare_results)
        fcn = self.start_opt
        for test in results.keys():
            failmsg = "Failure in "+test+"("+tstr(ttraj[test])+"): "
            self._single_test(fcn, ttraj[test], results[test], failmsg)

    def test_optional_start_can_append(self):
        bare_results = {'in' : True,
                        'out' : True,
                        'in_out' : True,
                        'out_in' : True,
                        'out_out_in' : True,
                        'in_out_in' : False,
                        'out_in_out' : True
                       }
        results = results_upper_lower(bare_results)
        fcn = self.start_opt.can_append
        for test in results.keys():
            failmsg = "Failure in "+test+"("+tstr(ttraj[test])+"): "
            self._single_test(fcn, ttraj[test], results[test], failmsg)

    def test_optional_start_can_preprend(self):
        bare_results = {'in' : True,
                        'out' : True,
                        'out_in_out' : True,
                        'out_out_in_out' : True,
                        'in_out' : True,
                        'out_in_out' : True,
                        'in_out_in_out' : False,
                        'out_in' : True,
                        'out_in_out_in' : False,
                        'in_out_in' : False
                       }
        results = results_upper_lower(bare_results)
        fcn = self.start_opt.can_prepend
        for test in results.keys():
            failmsg = "Failure in "+test+"("+tstr(ttraj[test])+"): "
            self._single_test(fcn, ttraj[test], results[test], failmsg)

    def test_optional_middle(self):
        bare_results = {'in_out_cross' : True,
                        'in_cross' :  True,
                        'in_out' : False,
                        'out_cross' : False,
                        'cross_in_cross_in' : False
                       }
        results = results_upper_lower(bare_results)
        fcn = self.mid_opt
        for test in results.keys():
            failmsg = "Failure in "+test+"("+tstr(ttraj[test])+"): "
            self._single_test(fcn, ttraj[test], results[test], failmsg)

    def test_optional_middle_can_append(self):
        bare_results = {'in' : True,
                        'out' : True,
                        'in_out' : True,
                        'out_in' : False,
                        'in_cross' : True,
                        'in_out_cross' : True,
                        'out_cross' : True,
                        'in_out_in' : False
                       }
        results = results_upper_lower(bare_results)
        fcn = self.mid_opt.can_append
        for test in results.keys():
            failmsg = "Failure in "+test+"("+tstr(ttraj[test])+"): "
            self._single_test(fcn, ttraj[test], results[test], failmsg)

    def test_optional_middle_can_preprend(self):
        bare_results = {'in' : True,
                        'out' : True,
                        'in_out' : True,
                        'out_in' : False,
                        'in_cross' : True,
                        'out_cross' : True,
                        'in_cross_in' : False
                       }
        results = results_upper_lower(bare_results)
        fcn = self.mid_opt.can_prepend
        for test in results.keys():
            failmsg = "Failure in "+test+"("+tstr(ttraj[test])+"): "
            self._single_test(fcn, ttraj[test], results[test], failmsg)

    def test_optional_end(self):
        bare_results = {'out_in' : True,
                        'out_in_out' : True,
                        'in_out' : False,
                        'out_out_in_out' : True,
                        'out_in_out_in' : False
                       }
        results = results_upper_lower(bare_results)
        fcn = self.end_opt
        for test in results.keys():
            failmsg = "Failure in "+test+"("+tstr(ttraj[test])+"): "
            self._single_test(fcn, ttraj[test], results[test], failmsg)

    def test_optional_end_can_append(self):
        bare_results = {'in' : True,
                        'out' : True,
                        'out_in' : True,
                        'in_out' : True,
                        'out_in_out' : True,
                        'in_out_in' : False,
                        'out_in_out_in' : False,
                        'in_in_out' : True
                       }
        results = results_upper_lower(bare_results)
        fcn = self.end_opt.can_append
        for test in results.keys():
            failmsg = "Failure in "+test+"("+tstr(ttraj[test])+"): "
            self._single_test(fcn, ttraj[test], results[test], failmsg)

    def test_optional_middle_can_prepend(self):
        bare_results = {'in' : True,
                        'out' : True,
                        'out_in' : True,
                        'in_out' : True,
                        'out_in_out' : True,
                        'in_out_in_out' : False,
                        'in_out_in' : False,
                        'out_in_out_in' : False
                       }
        results = results_upper_lower(bare_results)
        fcn = self.end_opt.can_prepend
        for test in results.keys():
            failmsg = "Failure in "+test+"("+tstr(ttraj[test])+"): "
            self._single_test(fcn, ttraj[test], results[test], failmsg)


    def test_optional_str(self):
        inX = InXEnsemble(vol1)
        opt_inX = OptionalEnsemble(inX)
        assert_equal(opt_inX.__str__(), "{"+inX.__str__()+"} (OPTIONAL)")

class testMinusInterfaceEnsemble(EnsembleTest):
    def setUp(self):
        # Mostly we use minus ensembles where the state matches the first
        # interface. We also test the case where that isn't in, in which
        # case there's an interstitial zone. (Only test it for nl=2 to keep
        # things easier.)
        self.minus_nl2 = MinusInterfaceEnsemble(
            state_vol=vol1,
            innermost_vol=vol1,
            n_l=2
        )
        self.minus_interstitial_nl2 = MinusInterfaceEnsemble(
            state_vol=vol1,
            innermost_vol=vol2,
            n_l=2
        )
        self.minus_nl3 = MinusInterfaceEnsemble(
            state_vol=vol1,
            innermost_vol=vol1,
            n_l=3
        )

    @raises(ValueError)
    def test_minus_nl1_fail(self):
        minus_nl1 = MinusInterfaceEnsemble(state_vol=vol1,
                                           innermost_vol=vol2,
                                           n_l=1)


    def test_minus_nl2_ensemble(self):
        non_default = [
            'in_cross_in_cross_in',
            'in_out_in_in_out_in',
            'in_out_in_out_in'
        ]
        self._test_everything(self.minus_nl2, non_default, False)

    def test_minus_nl2_can_append(self):
        non_default = [
            'in_cross_in_cross_in',
            'in_out_in_in_out_in',
            'in_out_in_out_in',
            'cross_in_cross_in',
            'in_in_cross_in',
            'in_in_out_in',
            'in_in_out_in_out',
            'in_in_out_out_in_in',
            'in_out_cross_out_in_out_in_out_cross_out_in',
            'out_in_cross_in',
            'out_in_in_in_out_in_out_in_in_in_out',
            'out_in_in_out_in',
            'out_in_out_in',
            'out_in_out_in_out',
            'out_in_out_out_in',
            'out_in_out_out_in_out',
            'in_hit_out_in_out',
            'out_hit_in_out_in'
        ]
        self._test_everything(self.minus_nl2.can_append, non_default, True)

    def test_minus_nl2_can_prepend(self):
        non_default = [
            'in_cross_in_cross',
            'in_cross_in_cross_in',
            'in_in_out_in_out',
            'in_in_out_out_in_in',
            'in_out_cross_out_in_out_in_out_cross_out_in',
            'in_out_in_in',
            'in_out_in_in_out',
            'in_out_in_in_out_in',
            'in_out_in_out',
            'in_out_in_out_in',
            'in_out_out_in_out',
            'out_in_in_in_out_in_out_in_in_in_out',
            'out_in_out_in_out',
            'out_in_out_out_in_out',
            'in_hit_out_in_out'
        ]
        self._test_everything(self.minus_nl2.can_prepend, non_default, True)

    def test_minus_interstitial_nl2_ensemble(self):
        non_default = [
            'in_cross_in_cross_in',
            'in_out_cross_out_in_out_in_out_cross_out_in',
        ]
        self._test_everything(self.minus_interstitial_nl2, non_default, False)

    def test_minus_interstitial_nl2_can_append(self):
        non_default = [
            'in_cross_in_cross_in',
            'in_out_cross_out_in_out_in_out_cross_out_in',
            'cross_in_cross_in',
            'in_in_cross_in',
            'out_in_cross_in'
        ]
        self._test_everything(self.minus_interstitial_nl2.can_append,
                              non_default, True)

    def test_minus_interstitial_nl2_can_prepend(self):
        non_default = [
            'in_cross_in_cross_in',
            'in_out_cross_out_in_out_in_out_cross_out_in',
            'in_cross_in_cross'
        ]
        self._test_everything(self.minus_interstitial_nl2.can_prepend,
                              non_default, True)

    def test_minus_nl3_ensemble(self):
        non_default = [
            'in_out_cross_out_in_out_in_out_cross_out_in',
        ]
        self._test_everything(self.minus_nl3, non_default, False)

    def test_minus_nl3_can_append(self):
        non_default = [
            'in_out_cross_out_in_out_in_out_cross_out_in',
            'out_in_in_in_out_in_out_in_in_in_out'
        ]
        self._test_everything(self.minus_nl3.can_append, non_default, True)

    def test_minus_nl3_can_prepend(self):
        non_default = [
            'in_out_cross_out_in_out_in_out_cross_out_in',
            'out_in_in_in_out_in_out_in_in_in_out'
        ]
        self._test_everything(self.minus_nl3.can_prepend, non_default, True)

