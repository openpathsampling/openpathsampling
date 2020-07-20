from __future__ import absolute_import
from builtins import zip
from builtins import map
from builtins import str
from builtins import range
from builtins import object
from nose.tools import (assert_equal, assert_not_equal, raises, assert_true,
                        assert_false)
from nose.plugins.skip import SkipTest
from .test_helpers import (CallIdentity, prepend_exception_message,
                          make_1d_traj, raises_with_message_like,
                          CalvinistDynamics)

import openpathsampling as paths
import openpathsampling.engines.openmm as peng
from openpathsampling.ensemble import *

import logging
logging.getLogger('openpathsampling.ensemble').setLevel(logging.DEBUG)
logging.getLogger('openpathsampling.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.netcdfplus').setLevel(logging.CRITICAL)
logger = logging.getLogger('openpathsampling.tests.testensemble')

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
            delta.append(adjustdict[char](random.randint(1, 3)))

        mydict[upperaddkey] = list(map(upper.__add__, delta))
        mydict[loweraddkey] = list(map(lower.__add__, delta))
        mydict[uppersubkey] = list(map(upper.__sub__, delta))
        mydict[lowersubkey] = list(map(lower.__sub__, delta))
    return mydict

def tstr(ttraj):
    return list(ttraj).__str__()

def results_upper_lower(adict):
    res_dict = {}
    for test in list(adict.keys()):
        res_dict['upper_'+test] = adict[test]
        res_dict['lower_'+test] = adict[test]
    return res_dict


def setup_module():
    ''' Setup for tests of classes in ensemble.py. '''
    #random.seed
    global lower, upper, op, vol1, vol2, vol3, ttraj, length0
    length0 = LengthEnsemble(0)
    lower = 0.1
    upper = 0.5
    op = paths.FunctionCV("Id", lambda snap : snap.coordinates[0][0])
    vol1 = paths.CVDefinedVolume(op, lower, upper).named('stateA')
    vol2 = paths.CVDefinedVolume(op, -0.1, 0.7).named('interface0')
    vol3 = paths.CVDefinedVolume(op, 2.0, 2.5).named('stateB')
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
                 "aabbaa", "abbb", "aaab"
                ]
    ttraj = build_trajdict(trajtypes, lower, upper)

    # make the tests from lists into trajectories
    for test in list(ttraj.keys()):
        ttraj[test] = make_1d_traj(coordinates=ttraj[test],
                                   velocities=[1.0]*len(ttraj[test]))

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
                to_append = 'out'
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

    def _test_everything(self, test_fcn, non_default=[], default=False):
        """
        Runs tests using *all* the trajectory test suite. This is the
        ultimate in test-running simplicity!!
        """
        results = {}
        for test in list(ttraj.keys()):
            results[test] = default
        nondef_dict = {}
        for test in non_default:
            if test in list(ttraj.keys()):
                results[test] = not default
            if "lower_"+test in list(ttraj.keys()):
                results["lower_"+test] = not default
            if "upper_"+test in list(ttraj.keys()):
                results["upper_"+test] = not default

        for test in list(results.keys()):
            logging.getLogger('openpathsampling.ensemble').debug(
                "Starting test for " + test + "("+str(ttraj[test])+")"
            )
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(test_fcn, ttraj[test], results[test], failmsg)


    def _run(self, results):
        """Actually run tests on the trajectory and the wrapped trajectory.

        Nearly all of the tests are just this simple. By creating custom error
        messages (using prepend_exception_message) we can wrap the many tests
        into loops instead of making tons of lines of code.
        """
        for test in list(results.keys()):
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.ensemble, ttraj[test], results[test], failmsg)

            wrapped = wrap_traj(ttraj[test], self.wrapstart, self.wrapend)
            lentt = len(ttraj[test])

            failmsg = "Failure in wrapped "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.ensemble, wrapped, results[test], failmsg)

            failmsg = "Failure in slice_ens "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.slice_ens, wrapped, results[test], failmsg)


class TestPartOutXEnsemble(EnsembleTest):
    def setup(self):
        self.leaveX = PartOutXEnsemble(vol1)

    def test_leaveX(self):
        """PartOutXEnsemble passes the trajectory test suite"""
        for test in list(ttraj.keys()):
            if "out" in in_out_parser(test):
                res = True
            else:
                res = False
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.leaveX, ttraj[test], res, failmsg)

    def test_invert(self):
        inverted = ~self.leaveX
        for test in list(ttraj.keys()):
            if "out" in in_out_parser(test):
                res = False
            else:
                res = True
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(inverted, ttraj[test], res, failmsg)

    def test_can_append(self):
        self._test_everything(self.leaveX.can_append, default=True)

    def test_can_prepend(self):
        self._test_everything(self.leaveX.can_prepend, default=True)

    def test_strict_can_append(self):
        self._test_everything(self.leaveX.strict_can_append, default=True)

    def test_strict_can_prepend(self):
        self._test_everything(self.leaveX.strict_can_prepend, default=True)

    def test_leaveX_0(self):
        """PartOutXEnsemble treatment of zero-length trajectory"""
        assert_equal(self.leaveX(paths.Trajectory([])), False)
        assert_equal(self.leaveX.can_append(paths.Trajectory([])), True)
        assert_equal(self.leaveX.can_prepend(paths.Trajectory([])), True)

    def test_leaveX_str(self):
        volstr = "{x|Id(x) in [0.1, 0.5]}"
        assert_equal(self.leaveX.__str__(), 
                     "exists t such that x[t] in (not "+volstr+")")

class TestAllInXEnsemble(EnsembleTest):
    def setup(self):
        self.inX = AllInXEnsemble(vol1)

    def test_inX(self):
        """AllInXEnsemble passes the trajectory test suite"""
        for test in list(ttraj.keys()):
            if "out" in in_out_parser(test):
                res = False
            else:
                res = True
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.inX, ttraj[test], res, failmsg)

    def test_can_append(self):
        for test in list(ttraj.keys()):
            if "out" in in_out_parser(test):
                res = False
            else:
                res = True
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.inX.can_append, ttraj[test], res, failmsg)

    def test_can_prepend(self):
        for test in list(ttraj.keys()):
            if "out" in in_out_parser(test):
                res = False
            else:
                res = True
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.inX.can_prepend, ttraj[test], res,
                              failmsg)

    def test_strict_can_append(self):
        for test in list(ttraj.keys()):
            if "out" in in_out_parser(test):
                res = False
            else:
                res = True
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.inX.strict_can_append, ttraj[test], res,
                              failmsg)

    def test_strict_can_prepend(self):
        for test in list(ttraj.keys()):
            if "out" in in_out_parser(test):
                res = False
            else:
                res = True
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.inX.strict_can_prepend, ttraj[test], res,
                              failmsg)

    def test_inX_0(self):
        """AllInXEnsemble treatment of zero-length trajectory"""
        assert_equal(self.inX(paths.Trajectory([])), False)
        assert_equal(self.inX.can_append(paths.Trajectory([])), True)
        assert_equal(self.inX.can_prepend(paths.Trajectory([])), True)

    def test_inX_str(self):
        volstr = "{x|Id(x) in [0.1, 0.5]}"
        assert_equal(self.inX.__str__(),
                     "x[t] in "+volstr+" for all t")

class TestAllOutXEnsemble(EnsembleTest):
    def setup(self):
        self.outX = AllOutXEnsemble(vol1)

    def test_outX(self):
        """AllOutXEnsemble passes the trajectory test suite"""
        for test in list(ttraj.keys()):
            if "in" in in_out_parser(test):
                res = False
            else:
                res = True
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.outX, ttraj[test], res, failmsg)

    def test_can_append(self):
        for test in list(ttraj.keys()):
            if "in" in in_out_parser(test):
                res = False
            else:
                res = True
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.outX.can_append, ttraj[test], res, failmsg)

    def test_can_prepend(self):
        for test in list(ttraj.keys()):
            if "in" in in_out_parser(test):
                res = False
            else:
                res = True
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.outX.can_prepend, ttraj[test], res, failmsg)

    def test_strict_can_append(self):
        for test in list(ttraj.keys()):
            if "in" in in_out_parser(test):
                res = False
            else:
                res = True
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.outX.strict_can_append, ttraj[test], res,
                              failmsg)

    def test_strict_can_prepend(self):
        for test in list(ttraj.keys()):
            if "in" in in_out_parser(test):
                res = False
            else:
                res = True
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.outX.strict_can_prepend, ttraj[test],
                              res, failmsg)

    def test_outX_0(self):
        """AllOutXEnsemble treatment of zero-length trajectory"""
        assert_equal(self.outX(paths.Trajectory([])), False)
        assert_equal(self.outX.can_append(paths.Trajectory([])), True)
        assert_equal(self.outX.can_prepend(paths.Trajectory([])), True)

    def test_outX_str(self):
        volstr = "{x|Id(x) in [0.1, 0.5]}"
        assert_equal(self.outX.__str__(),
                     "x[t] in (not "+volstr+") for all t")

class TestPartInXEnsemble(EnsembleTest):
    def setup(self):
        self.hitX = PartInXEnsemble(vol1)

    def test_hitX(self):
        """PartInXEnsemble passes the trajectory test suite"""
        for test in list(ttraj.keys()):
            if "in" in in_out_parser(test):
                res = True
            else:
                res = False
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.hitX, ttraj[test], res, failmsg)

    def test_can_append(self):
        self._test_everything(self.hitX.can_append, default=True)

    def test_can_prepend(self):
        self._test_everything(self.hitX.can_prepend, default=True)

    def test_strict_can_append(self):
        self._test_everything(self.hitX.strict_can_append, default=True)

    def test_strict_can_prepend(self):
        self._test_everything(self.hitX.strict_can_prepend, default=True)

    def test_hitX_0(self):
        """PartInXEnsemble treatment of zero-length trajectory"""
        assert_equal(self.hitX(paths.Trajectory([])), False)
        assert_equal(self.hitX.can_append(paths.Trajectory([])), True)
        assert_equal(self.hitX.can_prepend(paths.Trajectory([])), True)

    def test_hitX_str(self):
        volstr = "{x|Id(x) in [0.1, 0.5]}"
        assert_equal(self.hitX.__str__(),
                     "exists t such that x[t] in "+volstr)


class TestSequentialEnsemble(EnsembleTest):
    def setup(self):
        self.inX = AllInXEnsemble(vol1)
        self.outX = AllOutXEnsemble(vol1)
        self.hitX = PartInXEnsemble(vol1)
        self.leaveX = PartOutXEnsemble(vol1)
        self.inInterface = AllInXEnsemble(vol2)
        self.leaveX0 = PartOutXEnsemble(vol2)
        self.inX0 = AllInXEnsemble(vol2)
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

    def test_seqens_order_combo(self):
        # regression test for #229
        import numpy as np
        op = paths.FunctionCV(name="x", f=lambda snap : snap.xyz[0][0])
        bigvol = paths.CVDefinedVolume(collectivevariable=op,
                                    lambda_min=-100.0, lambda_max=100.0)

        traj = paths.Trajectory([
            paths.engines.toy.Snapshot(
                coordinates=np.array([[-0.5, 0.0]]),
                velocities=np.array([[0.0,0.0]])
            )
        ])

        vol_ens = paths.AllInXEnsemble(bigvol)
        len_ens = paths.LengthEnsemble(5)

        combo1 = vol_ens & len_ens
        combo2 = len_ens & vol_ens

        seq1 = SequentialEnsemble([combo1])
        seq2 = SequentialEnsemble([combo2])
        logger.debug("Checking combo1")
        assert_equal(combo1.can_append(traj), True)
        logger.debug("Checking combo2")
        assert_equal(combo2.can_append(traj), True)
        logger.debug("Checking seq1")
        assert_equal(seq1.can_append(traj), True)
        logger.debug("Checking seq2")
        assert_equal(seq2.can_append(traj), True)


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
        for test in list(results.keys()):
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.pseudo_tis.can_append, 
                                ttraj[test], results[test], failmsg)

    def test_strict_can_append_tis(self):
        results = {
            'upper_in_out' : True,
            'lower_in_out' : True,
            'upper_in_out_in' : False,
            'lower_in_out_in' : False,
            'upper_in' : True,
            'lower_in' : True,
            'upper_in_in_in' : False,
            'lower_in_in_in' : False,
            'upper_out_out_out' : False,
            'lower_out_out_out' : False,
            'upper_out_in' : False,
            'lower_out_in' : False,
            'upper_out' : False,
            'lower_out' : False,
            'upper_in_out_in_in' : False,
            'lower_in_out_in_in' : False,
            'upper_in_out_in_out_in' : False,
            'lower_in_out_in_out_in' : False
        }   
        for test in list(results.keys()):
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.pseudo_tis.strict_can_append, 
                                ttraj[test], results[test], failmsg)

    def test_can_append_pseudominus(self):
        """SequentialEnsemble as Pseudo-MinusEnsemble knows when it can append"""
        results = {
            'upper_in_out' : True,
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
        for test in list(results.keys()):
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.pseudo_minus.can_append, 
                                ttraj[test], results[test], failmsg)

    def test_strict_can_append_pseudominus(self):
        results = {
            'upper_in_out' : True,
            'lower_in_out' : True,
            'upper_in_out_in' : True,
            'lower_in_out_in' : True,
            'upper_in' : True,
            'lower_in' : True,
            'upper_in_in_in' : False,
            'lower_in_in_in' : False,
            'upper_out_out_out' : False,
            'lower_out_out_out' : False,
            'upper_out_in' : False,
            'lower_out_in' : False,
            'upper_out' : False,
            'lower_out' : False,
            'upper_in_out_in_in' : True,
            'lower_in_out_in_in' : True,
            'upper_in_out_in_out_in' : False,
            'lower_in_out_in_out_in' : False,
            'upper_in_out_in_in_out' : True,
            'lower_in_out_in_in_out' : True,
            'upper_out_in_out' : False,
            'lower_out_in_out' : False,
            'upper_out_in_in_out' : False,
            'lower_out_in_in_out' : False,
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
        for test in list(results.keys()):
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.pseudo_minus.strict_can_append, 
                                ttraj[test], results[test], failmsg)

    def test_can_prepend_pseudo_tis(self):
        """SequentialEnsemble as Pseudo-TISEnsemble knows when it can prepend"""
        results =   {
            'upper_in_out' : False,
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
        for test in list(results.keys()):
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.pseudo_tis.can_prepend, 
                                ttraj[test], results[test], failmsg)

    def test_strict_can_prepend_pseudo_tis(self):
        results =   {
            'upper_in_out' : False,
            'lower_in_out' : False,
            'upper_in_out_in' : False,
            'lower_in_out_in' : False,
            'upper_in' : True,
            'lower_in' : True,
            'upper_in_in_in' : False,
            'lower_in_in_in' : False,
            'upper_out_out_out' : False,
            'lower_out_out_out' : False,
            'upper_out_in' : True,
            'lower_out_in' : True,
            'upper_out' : False,
            'lower_out' : False,
            'upper_in_out_in_in' : False,
            'lower_in_out_in_in' : False,
            'upper_in_out_in_out_in' : False,
            'lower_in_out_in_out_in' : False
        }   
        for test in list(results.keys()):
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.pseudo_tis.strict_can_prepend, 
                                ttraj[test], results[test], failmsg)

    def test_can_prepend_pseudo_minus(self):
        results =   {
            'upper_in_out' : True,
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
        for test in list(results.keys()):
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.pseudo_minus.can_prepend, 
                                ttraj[test], results[test], failmsg)

    def test_strict_can_prepend_pseudo_minus(self):
        results =   {
            'upper_in_out' : False,
            'lower_in_out' : False,
            'upper_in_out_in' : True,
            'lower_in_out_in' : True,
            'upper_in' : True,
            'lower_in' : True,
            'upper_in_in_in' : False,
            'lower_in_in_in' : False,
            'upper_out_out_out' : False,
            'lower_out_out_out' : False,
            'upper_out_in' : True,
            'lower_out_in' : True,
            'upper_out' : False,
            'lower_out' : False,

            'upper_in_out_in_in' : False,
            'lower_in_out_in_in' : False,
            'upper_in_out_in_out_in' : False,
            'lower_in_out_in_out_in' : False,
            'upper_in_out_in_in_out' : False,
            'lower_in_out_in_in_out' : False,
            'upper_out_in_out' : False,
            'lower_out_in_out' : False,
            'upper_out_in_in_out' : False,
            'lower_out_in_in_out' : False,
            'upper_out_in_out_in': True,
            'lower_out_in_out_in': True,
            'upper_out_in_in_out_in' : True,
            'lower_out_in_in_out_in' : True
        }   
        for test in list(results.keys()):
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.pseudo_minus.strict_can_prepend, 
                                ttraj[test], results[test], failmsg)


    
    def test_sequential_transition_frames(self):
        """SequentialEnsemble identifies transitions frames correctly"""
        ensemble = SequentialEnsemble([self.inX, self.outX])
        results = {'upper_in_in_in' : [3],
                   'upper_out_out_out' : [],
                   'upper_in_out_in' : [1,2],
                   'upper_in_out' : [1,2]
                  }
        for test in list(results.keys()):
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
        for test in list(results.keys()):
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(ensemble, 
                                ttraj[test], results[test], failmsg)



    def test_sequential_in_out(self):
        """SequentialEnsembles based on In/AllOutXEnsemble"""
        # idea: for each ttraj, use the key name to define in/out behavior,
        # dynamically construct a SequentialEnsemble
        ens_dict = {'in' : self.inX, 'out' : self.outX }
        for test in list(ttraj.keys()):
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
        for test in list(ttraj.keys()):
            results[test] = False
        results['upper_in_out_in'] = True
        results['lower_in_out_in'] = True
        results['upper_in_out_out_in'] = True
        results['lower_in_out_out_in'] = True
        results['lower_in_out_cross_out_in'] = True
        results['upper_in_out_cross_out_in'] = True
        results['upper_in_cross_in'] = True
        results['lower_in_cross_in'] = True
        results['upper_in_hit_in'] = True
        for test in list(results.keys()):
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.pseudo_tis, ttraj[test], results[test],
                              failmsg)

    def test_sequential_pseudo_minus(self):
        """SequentialEnsemble as Pseudo-MinusEnsemble identifies paths"""
        results = {}
        for test in list(ttraj.keys()):
            results[test] = False
        results['upper_in_out_in_out_in'] = True
        results['lower_in_out_in_out_in'] = True
        results['upper_in_out_in_in_out_in'] = True
        results['lower_in_out_in_in_out_in'] = True
        results['upper_in_cross_in_cross_in'] = True
        results['lower_in_cross_in_cross_in'] = True
        for test in list(results.keys()):
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.pseudo_minus, ttraj[test], results[test],
                              failmsg)

    def test_sequential_tis(self):
        """SequentialEnsemble as TISEnsemble identifies paths"""
        results = {}
        for test in list(ttraj.keys()):
            results[test] = False
        results['upper_in_out_cross_out_in'] = True
        results['lower_in_out_cross_out_in'] = True
        results['upper_in_cross_in'] = True
        results['lower_in_cross_in'] = True
        for test in list(results.keys()):
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.tis, ttraj[test], results[test], failmsg)

    def test_sequential_generate_first_tis(self):
        """SequentialEnsemble to generate the first TIS path"""
        ensemble = SequentialEnsemble([
            self.outX | length0,
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
        logging.getLogger('openpathsampling.ensemble').info("Starting tests....")
        for test in list(match_results.keys()):
            failmsg = "Match failure in "+test+"("+str(ttraj[test])+"): "
            logging.getLogger('openpathsampling.ensemble').info(
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
        for test in list(append_results.keys()):
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


class TestSequentialEnsembleCombination(EnsembleTest):
    # testing EnsembleCombinations of SequentialEnsembles -- this is mainly
    # useful to making sure that the ensemble combination of strict_can_*
    # works correctly, since this is where strict and normal have a
    # distinction
    def setup(self):
        self.ens1 = SequentialEnsemble([
            AllInXEnsemble(vol1) & LengthEnsemble(1),
            AllOutXEnsemble(vol1) & PartOutXEnsemble(vol2),
            AllInXEnsemble(vol1) & LengthEnsemble(1)
        ])
        self.ens2 = SequentialEnsemble([
            AllInXEnsemble(vol1) & LengthEnsemble(1),
            LengthEnsemble(3),
            AllInXEnsemble(vol1) & LengthEnsemble(1)
        ])
        self.combo_and = self.ens1 & self.ens2
        self.combo_or = self.ens1 | self.ens2

    def test_call(self):
        ens1_passes = [
            'in_cross_in',
            'in_out_cross_in',
            'in_out_cross_out_in'
        ]
        self._test_everything(self.ens1, ens1_passes, False)
        ens2_passes = [
            'in_out_cross_out_in',
            'in_out_in_out_in',
            'in_cross_in_cross_in'
        ]
        self._test_everything(self.ens2, ens2_passes, False)

        or_passes = list(set(ens1_passes + ens2_passes))
        self._test_everything(self.combo_or, or_passes, False)

        and_passes = list(set(ens1_passes) & set(ens2_passes))
        self._test_everything(self.combo_and, and_passes, False)

    def test_can_append(self):
        ens1_true = [
            'hit',
            'in',
            'in_cross',
            'in_out',
            'in_out_cross',
            'in_out_out_out',
            'out',
            'out_cross',
            'out_out',
            'out_out_out',
            'upper_in_hit_out',
            'upper_out_hit_out'
        ]
        self._test_everything(self.ens1.can_append, ens1_true, False)
        ens2_true = [
            'hit',
            'in',
            'in_cross',
            'in_cross_in',
            'in_cross_in_cross',
            'in_hit_in',
            'in_in',
            'in_in_cross_in',
            'in_in_in',
            'in_in_in_out',
            'in_in_out',
            'in_in_out_in',
            'in_out',
            'in_out_cross',
            'in_out_in',
            'in_out_in_in',
            'in_out_in_out',
            'in_out_out_in',
            'in_out_out_out',
            'out',
            'out_cross',
            'out_hit_in',
            'out_in',
            'out_in_in',
            'out_in_out',
            'out_out',
            'out_out_in',
            'out_out_out',
            'upper_in_hit_out',
            'lower_in_hit_out',
            'upper_out_hit_out',
            'lower_out_hit_out'
        ]
        self._test_everything(self.ens2.can_append, ens2_true, False)

        or_true = list(set(ens1_true + ens2_true))
        self._test_everything(self.combo_or.can_append, or_true, False)

        and_true = list(set(ens1_true) & set(ens2_true))
        self._test_everything(self.combo_and.can_append, and_true, False)

    def test_can_prepend(self):
        ens1_true = [
            'hit',
            'in',
            'out',
            'out_cross',
            'out_in',
            'out_out',
            'out_out_in',
            'out_out_out',
            'out_out_out_in',
            'upper_out_hit_in',
            'upper_out_hit_out'
        ]
        self._test_everything(self.ens1.can_prepend, ens1_true, False)
        ens2_true = [
            'cross_in_cross_in',
            'hit',
            'in',
            'in_cross',
            'in_cross_in',
            'in_hit_in',
            'in_hit_out',
            'in_in',
            'in_in_cross_in',
            'in_in_in',
            'in_in_out',
            'in_in_out_in',
            'in_out',
            'in_out_cross',
            'in_out_in',
            'in_out_in_in',
            'in_out_out_in',
            'out',
            'out_cross',
            'upper_out_hit_in',
            'lower_out_hit_in',
            'upper_out_hit_out',
            'lower_out_hit_out',
            'out_in',
            'out_in_cross_in',
            'out_in_in',
            'out_in_in_in',
            'out_in_out',
            'out_in_out_in',
            'out_out',
            'out_out_in',
            'out_out_out',
            'out_out_out_in'
        ]
        self._test_everything(self.ens2.can_prepend, ens2_true, False)

        or_true = list(set(ens1_true + ens2_true))
        self._test_everything(self.combo_or.can_prepend, or_true, False)

        and_true = list(set(ens1_true) & set(ens2_true))
        self._test_everything(self.combo_and.can_prepend, and_true, False)

    def test_strict_can_append(self):
        ens1_true = [
            'lower_hit',
            'in',
            'in_cross',
            'in_out',
            'in_out_cross',
            'in_out_out_out',
            'upper_in_hit_out',
        ]
        self._test_everything(self.ens1.strict_can_append, ens1_true, False)
        ens2_true = [
            'lower_hit',
            'in',
            'in_cross',
            'in_cross_in',
            'in_cross_in_cross',
            'in_hit_in',
            'upper_in_hit_out',
            'lower_in_hit_out',
            'in_in',
            'in_in_cross_in',
            'in_in_in',
            'in_in_in_out',
            'in_in_out',
            'in_in_out_in',
            'in_out',
            'in_out_cross',
            'in_out_in',
            'in_out_in_in',
            'in_out_in_out',
            'in_out_out_in',
            'in_out_out_out',
        ]
        self._test_everything(self.ens2.strict_can_append, ens2_true, False)

        or_true = list(set(ens1_true + ens2_true))
        self._test_everything(self.combo_or.strict_can_append, or_true, False)

        and_true = list(set(ens1_true) & set(ens2_true))
        self._test_everything(self.combo_and.strict_can_append, and_true, False)

    def test_strict_can_prepend(self):
        ens1_true = [
            'lower_hit',
            'in',
            'out_in',
            'out_out_in',
            'out_out_out_in',
            'upper_out_hit_in',
        ]
        self._test_everything(self.ens1.strict_can_prepend, ens1_true, False)
        ens2_true = [
            'cross_in_cross_in',
            'lower_hit',
            'in',
            'in_cross_in',
            'in_hit_in',
            'in_in',
            'in_in_cross_in',
            'in_in_in',
            'in_in_out_in',
            'in_out_in',
            'in_out_in_in',
            'in_out_out_in',
            'upper_out_hit_in',
            'lower_out_hit_in',
            'out_in',
            'out_in_cross_in',
            'out_in_in',
            'out_in_in_in',
            'out_in_out_in',
            'out_out_in',
            'out_out_out_in'
        ]
        self._test_everything(self.ens2.strict_can_prepend, ens2_true, False)

        or_true = list(set(ens1_true + ens2_true))
        self._test_everything(self.combo_or.strict_can_prepend, or_true, False)

        and_true = list(set(ens1_true) & set(ens2_true))
        self._test_everything(self.combo_and.strict_can_prepend, and_true, False)

class TestTISEnsemble(EnsembleTest):
    def setup(self):
        self.tis = TISEnsemble(vol1, vol3, vol2, op)
        self.traj = ttraj['upper_in_out_cross_out_in']
        self.minl = min(op(self.traj))
        self.maxl = max(op(self.traj))

    def test_tis_trajectory_summary(self):
        summ = self.tis.trajectory_summary(self.traj)
        assert_equal(summ['initial_state'], 0)
        assert_equal(summ['final_state'], 0)
        assert_equal(summ['max_lambda'], self.maxl)
        assert_equal(summ['min_lambda'], self.minl)

    def test_tis_trajectory_summary_str(self):
        mystr = self.tis.trajectory_summary_str(self.traj)
        teststr = ("initial_state=stateA final_state=stateA min_lambda=" +
                   str(self.minl) + " max_lambda=" + str(self.maxl) + " ")
        assert_equal(mystr, teststr)

    def test_no_frame_after_interface(self):
        traj_3 = make_1d_traj([0.2, 0.6,  2.1])
        assert_equal(self.tis(traj_3), True)
        traj_2 = make_1d_traj([0.2, 2.1])
        assert_equal(self.tis(traj_2), True)

    def test_tis_ensemble_candidate(self):
        tis = TISEnsemble(vol1, vol3, vol2, op, lambda_i=0.7)
        test_f = lambda t: tis(t, candidate=True)
        results = {}
        upper_keys = [k for k in list(ttraj.keys()) if k[:6] == "upper_"]
        for test in upper_keys:
            results[test] = False
        # results where we SHOULD get True
        results['upper_in_out_cross_out_in'] = True
        results['upper_in_cross_in'] = True
        # results where the input isn't actually a candidate trajectory, so
        # it accepts the path even though it shouldn't
        results['upper_in_cross_in_cross_in'] = True
        results['upper_in_out_cross_out_in_out_in_out_cross_out_in'] = True
        results['upper_in_in_cross_in'] = True

        for test in list(results.keys()):
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(test_f, ttraj[test], results[test], failmsg)

    def test_tis_ensemble_candidate_cv_max(self):
        cv_max_func = lambda t, cv_: max(cv_(t))
        cv_max = paths.netcdfplus.FunctionPseudoAttribute(
            name="max " + op.name,
            key_class=paths.Trajectory,
            f=lambda t: cv_max_func(t, cv_=op)
        )
        tis = TISEnsemble(vol1, vol3, vol2, cv_max=cv_max, lambda_i=0.7)

        test_f = lambda t: tis(t, candidate=True)
        results = {}
        upper_keys = [k for k in list(ttraj.keys()) if k[:6] == "upper_"]
        for test in upper_keys:
            results[test] = False
        # results where we SHOULD get True
        results['upper_in_out_cross_out_in'] = True
        results['upper_in_cross_in'] = True
        # results where the input isn't actually a candidate trajectory, so
        # it accepts the path even though it shouldn't
        results['upper_in_cross_in_cross_in'] = True
        results['upper_in_out_cross_out_in_out_in_out_cross_out_in'] = True
        results['upper_in_in_cross_in'] = True

        # first test should cache the results
        for test in list(results.keys()):
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(test_f, ttraj[test], results[test], failmsg)

        # second test should check the cache
        for test in list(results.keys()):
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(test_f, ttraj[test], results[test], failmsg)


class EnsembleCacheTest(EnsembleTest):
    def _was_cache_reset(self, cache):
        return cache.contents == { }

class TestEnsembleCache(EnsembleCacheTest):
    def setup(self):
        self.fwd = EnsembleCache(direction=+1)
        self.rev = EnsembleCache(direction=-1)
        self.traj = ttraj['lower_in_out_in_in_out_in']

    def test_initially_reset(self):
        assert_equal(self._was_cache_reset(self.fwd), True)
        assert_equal(self._was_cache_reset(self.rev), True)

    def test_change_trajectory(self):
        traj2 = ttraj['lower_in_out_in']
        # tests for forward
        self.fwd.contents = { 'test' : 'object' }
        assert_equal(self._was_cache_reset(self.fwd), False)
        self.fwd.check(self.traj)
        assert_equal(self._was_cache_reset(self.fwd), True)
        self.fwd.contents['ens_num'] = 1
        assert_equal(self._was_cache_reset(self.fwd), False)
        self.fwd.check(traj2)
        assert_equal(self._was_cache_reset(self.fwd), True)
        # tests for backward
        self.rev.contents = { 'test' : 'object' }
        assert_equal(self._was_cache_reset(self.rev), False)
        self.rev.check(self.traj)
        assert_equal(self._was_cache_reset(self.rev), True)
        self.rev.contents['ens_num'] = 1
        assert_equal(self._was_cache_reset(self.rev), False)
        self.rev.check(traj2)
        assert_equal(self._was_cache_reset(self.rev), True)

    def test_trajectory_by_frame(self):
        # tests for forward
        self.fwd.check(self.traj[0:1])
        assert_equal(self._was_cache_reset(self.fwd), True)
        self.fwd.contents = { 'test' : 'object' }
        assert_equal(self._was_cache_reset(self.fwd), False)
        self.fwd.check(self.traj[0:2])
        assert_equal(self._was_cache_reset(self.fwd), False)
        # tests for backward
        self.rev.check(self.traj[-1:])
        assert_equal(self._was_cache_reset(self.rev), True)
        self.rev.contents = { 'test' : 'object' }
        assert_equal(self._was_cache_reset(self.rev), False)
        self.rev.check(self.traj[-2:])
        assert_equal(self._was_cache_reset(self.rev), False)

    def test_same_traj_twice_no_reset(self):
        # tests for forward
        self.fwd.check(self.traj)
        assert_equal(self._was_cache_reset(self.fwd), True)
        self.fwd.contents = { 'test' : 'object' }
        self.fwd.check(self.traj)
        assert_equal(self._was_cache_reset(self.fwd), False)
        # tests for backward
        self.rev.check(self.traj)
        assert_equal(self._was_cache_reset(self.rev), True)
        self.rev.contents = { 'test' : 'object' }
        self.rev.check(self.traj)
        assert_equal(self._was_cache_reset(self.rev), False)


    def test_trajectory_skips_frame(self):
        # tests for forward
        self.fwd.check(self.traj[0:1])
        assert_equal(self._was_cache_reset(self.fwd), True)
        self.fwd.contents = { 'test' : 'object' }
        assert_equal(self._was_cache_reset(self.fwd), False)
        self.fwd.check(self.traj[0:3])
        assert_equal(self._was_cache_reset(self.fwd), True)
        # tests for backward
        self.rev.check(self.traj[-1:])
        assert_equal(self._was_cache_reset(self.rev), True)
        self.rev.contents = { 'test' : 'object' }
        assert_equal(self._was_cache_reset(self.rev), False)
        self.rev.check(self.traj[-3:])
        assert_equal(self._was_cache_reset(self.rev), True)

    def test_trajectory_middle_frame_changes(self):
        # tests for forward
        self.fwd.check(self.traj[0:2])
        assert_equal(self._was_cache_reset(self.fwd), True)
        self.fwd.contents = { 'test' : 'object' }
        assert_equal(self._was_cache_reset(self.fwd), False)
        new_traj = self.traj[0:1] + self.traj[3:5]
        self.fwd.check(new_traj)
        assert_equal(self._was_cache_reset(self.fwd), True)
        # tests for backward
        self.rev.check(self.traj[0:2])
        assert_equal(self._was_cache_reset(self.rev), True)
        self.rev.contents = { 'test' : 'object' }
        assert_equal(self._was_cache_reset(self.rev), False)
        new_traj = self.traj[-4:-2] + self.traj[-1:]
        self.rev.check(new_traj)
        assert_equal(self._was_cache_reset(self.rev), True)


class TestSequentialEnsembleCache(EnsembleCacheTest):
    def setup(self):
        self.inX = AllInXEnsemble(vol1)
        self.outX = AllOutXEnsemble(vol1)
        self.length1 = LengthEnsemble(1)
        self.pseudo_minus = SequentialEnsemble([
            self.inX & self.length1,
            self.outX,
            self.inX,
            self.outX,
            self.inX & self.length1
        ])
        self.traj = ttraj['lower_in_out_in_in_out_in']

    def test_all_in_as_seq_can_append(self):
        ens = SequentialEnsemble([AllInXEnsemble(vol1 | vol2 | vol3)])
        cache = ens._cache_can_append
        traj = ttraj['upper_in_in_out_out_in_in']
        assert_equal(ens.can_append(traj[0:1]), True)
        assert_equal(ens.can_append(traj[0:2]), True)
        assert_equal(ens.can_append(traj[0:3]), True)
        assert_equal(ens.can_append(traj[0:4]), True)
        assert_equal(ens.can_append(traj[0:5]), True)
        assert_equal(ens.can_append(traj[0:6]), True)

    def test_sequential_caching_can_append(self):
        cache = self.pseudo_minus._cache_can_append
        assert_equal(self.pseudo_minus.can_append(self.traj[0:1]), True)
        assert_equal(cache.contents['ens_num'], 1)
        assert_equal(cache.contents['ens_from'], 0)
        assert_equal(cache.contents['subtraj_from'], 1)
        logging.getLogger('openpathsampling.ensemble').debug("Starting [0:2]")
        assert_equal(self.pseudo_minus.can_append(self.traj[0:2]), True)
        assert_equal(cache.contents['ens_num'], 1)
        assert_equal(cache.contents['ens_from'], 0)
        assert_equal(cache.contents['subtraj_from'], 1)
        logging.getLogger('openpathsampling.ensemble').debug("Starting [0:3]")
        assert_equal(self.pseudo_minus.can_append(self.traj[0:3]), True)
        assert_equal(cache.contents['ens_num'], 2)
        assert_equal(cache.contents['ens_from'], 0)
        assert_equal(cache.contents['subtraj_from'], 2)
        logging.getLogger('openpathsampling.ensemble').debug("Starting [0:4]")
        assert_equal(self.pseudo_minus.can_append(self.traj[0:4]), True)
        assert_equal(cache.contents['ens_num'], 2)
        assert_equal(cache.contents['ens_from'], 0)
        assert_equal(cache.contents['subtraj_from'], 2)
        logging.getLogger('openpathsampling.ensemble').debug("Starting [0:5]")
        assert_equal(self.pseudo_minus.can_append(self.traj[0:5]), True)
        assert_equal(cache.contents['ens_num'], 3)
        assert_equal(cache.contents['ens_from'], 0)
        assert_equal(cache.contents['subtraj_from'], 4)
        logging.getLogger('openpathsampling.ensemble').debug("Starting [0:6]")
        assert_equal(self.pseudo_minus.can_append(self.traj[0:6]), False)
        assert_equal(cache.contents['ens_num'], 4)
        assert_equal(cache.contents['ens_from'], 0)
        assert_equal(cache.contents['subtraj_from'], 5)

    def test_sequential_caching_resets(self):
        #cache = self.pseudo_minus._cache_can_append
        assert_equal(self.pseudo_minus.can_append(self.traj[2:3]), True)
        assert_equal(self.pseudo_minus(self.traj[2:3]), False)
        #assert_equal(self._was_cache_reset(cache), True)
        assert_equal(self.pseudo_minus.can_append(self.traj[2:4]), True)
        assert_equal(self.pseudo_minus(self.traj[2:4]), False)
        #assert_equal(self._was_cache_reset(cache), True)
        for i in range(4, len(self.traj)-1):
            assert_equal(self.pseudo_minus.can_append(self.traj[2:i+1]), True)
            assert_equal(self.pseudo_minus(self.traj[2:i+1]), False)
            #assert_equal(self._was_cache_reset(cache), False)
        assert_equal(self.pseudo_minus.can_append(self.traj[2:]), False)
        assert_equal(self.pseudo_minus(self.traj[2:]), False)
        #assert_equal(self._was_cache_reset(cache), False)
        # TODO: same story backward
        raise SkipTest

    def test_sequential_caching_call(self):
        raise SkipTest

    def test_sequential_caching_can_prepend(self):
        cache = self.pseudo_minus._cache_can_prepend
        assert_equal(self.pseudo_minus.can_prepend(self.traj[5:6]), True)
        assert_equal(cache.contents['ens_num'], 3)
        assert_equal(cache.contents['ens_from'], 4)
        assert_equal(cache.contents['subtraj_from'], -1)
        assert_equal(self.pseudo_minus.can_prepend(self.traj[4:6]), True)
        assert_equal(cache.contents['ens_num'], 3)
        assert_equal(cache.contents['ens_from'], 4)
        assert_equal(cache.contents['subtraj_from'], -1)
        assert_equal(self.pseudo_minus.can_prepend(self.traj[3:6]), True)
        assert_equal(cache.contents['ens_num'], 2)
        assert_equal(cache.contents['ens_from'], 4)
        assert_equal(cache.contents['subtraj_from'], -2)
        assert_equal(self.pseudo_minus.can_prepend(self.traj[2:6]), True)
        assert_equal(cache.contents['ens_num'], 2)
        assert_equal(cache.contents['ens_from'], 4)
        assert_equal(cache.contents['subtraj_from'], -2)
        assert_equal(self.pseudo_minus.can_prepend(self.traj[1:6]), True)
        assert_equal(cache.contents['ens_num'], 1)
        assert_equal(cache.contents['ens_from'], 4)
        assert_equal(cache.contents['subtraj_from'], -4)
        assert_equal(self.pseudo_minus.can_prepend(self.traj[0:6]), False)
        assert_equal(cache.contents['ens_num'], 0)
        assert_equal(cache.contents['ens_from'], 4)
        assert_equal(cache.contents['subtraj_from'], -5)




class TestSlicedTrajectoryEnsemble(EnsembleTest):
    def test_sliced_ensemble_init(self):
        init_as_int = SlicedTrajectoryEnsemble(AllInXEnsemble(vol1), 3)
        init_as_slice = SlicedTrajectoryEnsemble(AllInXEnsemble(vol1),
                                                 slice(3, 4))
        assert_equal(init_as_int, init_as_slice)
        assert_equal(init_as_slice.region, init_as_int.region)

    def test_sliced_as_TISEnsemble(self):
        '''SlicedTrajectory and Sequential give same TIS results'''
        sliced_tis = (
            SlicedTrajectoryEnsemble(AllInXEnsemble(vol1), 0) &
            SlicedTrajectoryEnsemble(AllOutXEnsemble(vol1 | vol3), slice(1,-1)) &
            SlicedTrajectoryEnsemble(PartOutXEnsemble(vol2), slice(1,-1)) &
            SlicedTrajectoryEnsemble(AllInXEnsemble(vol1 | vol3), -1)
        )
        sequential_tis = SequentialEnsemble([
            AllInXEnsemble(vol1) & LengthEnsemble(1),
            AllOutXEnsemble(vol1 | vol3) & PartOutXEnsemble(vol2),
            AllInXEnsemble(vol1 | vol3) & LengthEnsemble(1)
        ])
        real_tis = paths.TISEnsemble(vol1, vol3, vol2)
        for test in list(ttraj.keys()):
            failmsg = "Failure in "+test+"("+tstr(ttraj[test])+"): "
            self._single_test(real_tis, ttraj[test],
                              sequential_tis(ttraj[test]), failmsg)
            self._single_test(sliced_tis, ttraj[test], 
                              sequential_tis(ttraj[test]), failmsg)

    def test_slice_outside_trajectory_range(self):
        ens = SlicedTrajectoryEnsemble(AllInXEnsemble(vol1), slice(5,9))
        test = 'upper_in'
        # the slice should return the empty trajectory, and therefore should
        # return false
        assert_equal(ens(ttraj[test]), False)

    def test_even_sliced_trajectory(self):
        even_slice = slice(None, None, 2)
        ens = SlicedTrajectoryEnsemble(AllInXEnsemble(vol1), even_slice)
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
        for test in list(results.keys()):
            failmsg = "Failure in "+test+"("+tstr(ttraj[test])+"): "
            self._single_test(ens, ttraj[test], results[test], failmsg)

    def test_sliced_sequential_global_whole(self):
        even_slice = slice(None, None, 2)
        ens = SlicedTrajectoryEnsemble(SequentialEnsemble([
            AllInXEnsemble(vol1),
            AllOutXEnsemble(vol1)
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
        for test in list(results.keys()):
            failmsg = "Failure in "+test+"("+tstr(ttraj[test])+"): "
            self._single_test(ens, ttraj[test], results[test], failmsg)

    def test_sliced_sequential_subtraj_member(self):
        even_slice = slice(None, None, 2)
        ens = SequentialEnsemble([
            AllInXEnsemble(vol1),
            SlicedTrajectoryEnsemble(AllOutXEnsemble(vol1), even_slice)
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
        for test in list(results.keys()):
            failmsg = "Failure in "+test+"("+tstr(ttraj[test])+"): "
            self._single_test(ens, ttraj[test], results[test], failmsg)

    def test_sliced_sequential_subtraj_middle(self):
        even_slice = slice(None, None, 2)
        ens = SequentialEnsemble([
            AllInXEnsemble(vol1),
            SlicedTrajectoryEnsemble(AllOutXEnsemble(vol1), even_slice),
            AllInXEnsemble(vol1) & LengthEnsemble(1)
        ])
        bare_results = {'in_in_out_out_in_in' : False
                       }
        results = results_upper_lower(bare_results)
        for test in list(results.keys()):
            failmsg = "Failure in "+test+"("+tstr(ttraj[test])+"): "
            self._single_test(ens, ttraj[test], results[test], failmsg)


    def test_sliced_str(self):
        even_slice = slice(None,None, 2)
        slice_1_10 = slice(1, 10)
        slice_1_end = slice(1,None)
        slice_no_ends = slice(1, -1)
        inX = AllInXEnsemble(vol1)
        inXstr = "x[t] in {x|Id(x) in [0.1, 0.5]} for all t"
        assert_equal(SlicedTrajectoryEnsemble(inX, even_slice).__str__(),
                     "("+inXstr+" in {:} every 2)")
        assert_equal(SlicedTrajectoryEnsemble(inX, slice_1_10).__str__(),
                     "("+inXstr+" in {1:10})")
        assert_equal(SlicedTrajectoryEnsemble(inX, slice_1_end).__str__(),
                     "("+inXstr+" in {1:})")
        assert_equal(SlicedTrajectoryEnsemble(inX, slice_no_ends).__str__(),
                     "("+inXstr+" in {1:-1})")

class TestOptionalEnsemble(EnsembleTest):
    def setup(self):
        self.start_opt = SequentialEnsemble([
            OptionalEnsemble(AllOutXEnsemble(vol1)),
            AllInXEnsemble(vol1),
            AllOutXEnsemble(vol1)
        ])
        self.end_opt = SequentialEnsemble([
            AllOutXEnsemble(vol1),
            AllInXEnsemble(vol1),
            OptionalEnsemble(AllOutXEnsemble(vol1))
        ])
        self.mid_opt = SequentialEnsemble([
            AllInXEnsemble(vol1),
            OptionalEnsemble(AllOutXEnsemble(vol1) & AllInXEnsemble(vol2)),
            AllOutXEnsemble(vol2)
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
        for test in list(results.keys()):
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
        for test in list(results.keys()):
            failmsg = "Failure in "+test+"("+tstr(ttraj[test])+"): "
            self._single_test(fcn, ttraj[test], results[test], failmsg)

    def test_optional_start_strict_can_append(self):
        bare_results = {'in' : True,
                        'out' : True,
                        'in_out' : True,
                        'out_in' : True,
                        'out_out_in' : True,
                        'in_out_in' : False,
                        'out_in_out' : True
                       }
        results = results_upper_lower(bare_results)
        fcn = self.start_opt.strict_can_append
        for test in list(results.keys()):
            failmsg = "Failure in "+test+"("+tstr(ttraj[test])+"): "
            self._single_test(fcn, ttraj[test], results[test], failmsg)

    def test_optional_start_can_prepend(self):
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
        for test in list(results.keys()):
            failmsg = "Failure in "+test+"("+tstr(ttraj[test])+"): "
            self._single_test(fcn, ttraj[test], results[test], failmsg)

    def test_optional_start_strict_can_prepend(self):
        bare_results = {'in' : False,
                        'out' : True,
                        'out_in_out' : True,
                        'out_out_in_out' : True,
                        'in_out' : True,
                        'out_in_out' : True,
                        'in_out_in_out' : False,
                        'out_in' : False,
                        'out_in_out_in' : False,
                        'in_out_in' : False
                       }
        results = results_upper_lower(bare_results)
        fcn = self.start_opt.strict_can_prepend
        for test in list(results.keys()):
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
        for test in list(results.keys()):
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
        for test in list(results.keys()):
            failmsg = "Failure in "+test+"("+tstr(ttraj[test])+"): "
            self._single_test(fcn, ttraj[test], results[test], failmsg)

    def test_optional_middle_strict_can_append(self):
        bare_results = {'in' : True,
                        'out' : False,
                        'in_out' : True,
                        'out_in' : False,
                        'in_cross' : True,
                        'in_out_cross' : True,
                        'out_cross' : False,
                        'in_out_in' : False
                       }
        results = results_upper_lower(bare_results)
        fcn = self.mid_opt.strict_can_append
        for test in list(results.keys()):
            failmsg = "Failure in "+test+"("+tstr(ttraj[test])+"): "
            self._single_test(fcn, ttraj[test], results[test], failmsg)

    def test_optional_middle_can_prepend(self):
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
        for test in list(results.keys()):
            failmsg = "Failure in "+test+"("+tstr(ttraj[test])+"): "
            self._single_test(fcn, ttraj[test], results[test], failmsg)

    def test_optional_middle_strict_can_prepend(self):
        bare_results = {'in' : False,
                        'out' : False,
                        'in_out' : False,
                        'out_in' : False,
                        'in_cross' : True,
                        'out_cross' : True,
                        'in_cross_in' : False
                       }
        results = results_upper_lower(bare_results)
        fcn = self.mid_opt.strict_can_prepend
        for test in list(results.keys()):
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
        for test in list(results.keys()):
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
        for test in list(results.keys()):
            failmsg = "Failure in "+test+"("+tstr(ttraj[test])+"): "
            self._single_test(fcn, ttraj[test], results[test], failmsg)

    def test_optional_end_strict_can_append(self):
        bare_results = {'in' : False,
                        'out' : True,
                        'out_in' : True,
                        'in_out' : False,
                        'out_in_out' : True,
                        'in_out_in' : False,
                        'out_in_out_in' : False,
                        'in_in_out' : False
                       }
        results = results_upper_lower(bare_results)
        fcn = self.end_opt.strict_can_append
        for test in list(results.keys()):
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
        for test in list(results.keys()):
            failmsg = "Failure in "+test+"("+tstr(ttraj[test])+"): "
            self._single_test(fcn, ttraj[test], results[test], failmsg)

    def test_optional_middle_strict_can_prepend(self):
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
        fcn = self.end_opt.strict_can_prepend
        for test in list(results.keys()):
            failmsg = "Failure in "+test+"("+tstr(ttraj[test])+"): "
            self._single_test(fcn, ttraj[test], results[test], failmsg)


    def test_optional_str(self):
        inX = AllInXEnsemble(vol1)
        opt_inX = OptionalEnsemble(inX)
        assert_equal(opt_inX.__str__(), "{"+inX.__str__()+"} (OPTIONAL)")

class TestPrefixTrajectoryEnsemble(EnsembleTest):
    def setup(self):
        self.inX = AllInXEnsemble(vol1)

    def test_bad_start_traj(self):
        traj = ttraj['upper_out_in_in_in']
        ens = PrefixTrajectoryEnsemble(
            SequentialEnsemble([self.inX]),
            traj[0:2]
        )
        assert_equal(ens.can_append(traj[0:3]), False)
        assert_equal(ens.strict_can_append(traj[0:3]), False)
        assert_equal(ens(traj[0:3]), False)

    def test_good_start_traj(self):
        traj = ttraj['upper_in_in_in']
        ens = PrefixTrajectoryEnsemble(
            SequentialEnsemble([self.inX]),
            traj[0:2]
        )
        assert_equal(ens.can_append(traj[2:3]), True)
        assert_equal(ens.strict_can_append(traj[2:3]), True)
        assert_equal(ens(traj[2:3]), True)

    @raises(RuntimeError)
    def test_can_prepend(self):
        traj = ttraj['upper_in_in_in']
        ens = PrefixTrajectoryEnsemble(
            SequentialEnsemble([self.inX]),
            traj[0:2]
        )
        ens.can_prepend(traj[2:3])

    @raises(RuntimeError)
    def test_strict_can_prepend(self):
        traj = ttraj['upper_in_in_in']
        ens = PrefixTrajectoryEnsemble(
            SequentialEnsemble([self.inX]),
            traj[0:2]
        )
        ens.strict_can_prepend(traj[2:3])

    def test_caching_in_fwdapp_seq(self):
        inX = AllInXEnsemble(vol1)
        outX = AllOutXEnsemble(vol1)
        length1 = LengthEnsemble(1)
        pseudo_minus = SequentialEnsemble([
            inX & length1,
            outX,
            inX,
            outX,
            inX & length1 
        ])
        traj = ttraj['upper_in_out_in_in_out_in']
        ens = PrefixTrajectoryEnsemble(pseudo_minus, traj[0:2])
        assert_equal(ens.can_append(traj[2:3]), True)
        assert_equal(ens._cached_trajectory, traj[0:3])
        assert_equal(ens._cache_can_append.trusted, False)

        assert_equal(ens.can_append(traj[2:4]), True)
        assert_equal(ens._cached_trajectory, traj[0:4])
        assert_equal(ens._cache_can_append.trusted, True)

        assert_equal(ens.can_append(traj[2:5]), True)
        assert_equal(ens._cached_trajectory, traj[0:5])
        assert_equal(ens._cache_can_append.trusted, True)

        assert_equal(ens.can_append(traj[2:6]), False)
        assert_equal(ens._cached_trajectory, traj[0:6])
        assert_equal(ens._cache_can_append.trusted, True)


class TestSuffixTrajectoryEnsemble(EnsembleTest):
    def setup(self):
        xval = paths.FunctionCV("x", lambda s : s.xyz[0][0])
        vol = paths.CVDefinedVolume(xval, 0.1, 0.5)
        self.inX = AllInXEnsemble(vol)
        self.outX = AllOutXEnsemble(vol)

    def test_bad_end_traj(self):
        traj = ttraj['upper_in_in_in_out']
        ens = SuffixTrajectoryEnsemble(
            SequentialEnsemble([self.inX]),
            traj[-2:]
        )
        assert_equal(ens.can_prepend(traj[-3:2]), False)
        assert_equal(ens.strict_can_prepend(traj[-3:2]), False)
        assert_equal(ens(traj[-3:2]), False)

    def test_good_end_traj(self):
        traj = ttraj['upper_out_in_in_in']
        ens = SuffixTrajectoryEnsemble(
            SequentialEnsemble([self.inX]),
            traj[-2:]
        )
        assert_equal(ens.can_prepend(traj[-3:-2]), True)
        assert_equal(ens.strict_can_prepend(traj[-3:-2]), True)
        assert_equal(ens(traj[-3:-2]), True)
        assert_equal(ens.can_prepend(traj[-4:-2]), False)
        assert_equal(ens.strict_can_prepend(traj[-4:-2]), False)
        assert_equal(ens(traj[-4:-2]), False)

    @raises(RuntimeError)
    def test_can_append(self):
        traj = ttraj['upper_out_in_in_in']
        ens = SuffixTrajectoryEnsemble(
            SequentialEnsemble([self.inX]),
            traj[-2:]
        )
        ens.can_append(traj)

    @raises(RuntimeError)
    def test_strict_can_append(self):
        traj = ttraj['upper_out_in_in_in']
        ens = SuffixTrajectoryEnsemble(
            SequentialEnsemble([self.inX]),
            traj[-2:]
        )
        ens.strict_can_append(traj)

    def test_caching_in_bkwdprep_seq(self):
        length1 = LengthEnsemble(1)
        pseudo_minus = SequentialEnsemble([
            self.inX & length1,
            self.outX,
            self.inX,
            self.outX,
            self.inX & length1 
        ])
        traj = ttraj['upper_in_out_in_in_out_in']

        # sanity checks before running the suffixed version
        assert_equal(pseudo_minus(traj), True)
        for i in range(-1, -6):
            assert_equal(pseudo_minus.can_prepend(traj[i:]), True)

        logger.debug("alltraj " + str([id(i) for i in traj]))
        ens = SuffixTrajectoryEnsemble(pseudo_minus, traj[-3:])
        assert_equal(len(ens._cached_trajectory), 3)

        assert_equal(ens.can_prepend(traj[-4:-3].reversed), True)
        assert_equal(len(ens._cached_trajectory), 4)
        assert_equal(ens._cache_can_prepend.trusted, False)

        assert_equal(ens.can_prepend(traj[-5:-3].reversed), True)
        assert_equal(len(ens._cached_trajectory), 5)
        assert_equal(ens._cache_can_prepend.trusted, True)

        assert_equal(ens.can_prepend(traj[-6:-3].reversed), False)
        assert_equal(len(ens._cached_trajectory), 6)
        assert_equal(ens._cache_can_prepend.trusted, True)

class TestMinusInterfaceEnsemble(EnsembleTest):
    def setup(self):
        # Mostly we use minus ensembles where the state matches the first
        # interface. We also test the case where that isn't in, in which
        # case there's an interstitial zone. (Only test it for nl=2 to keep
        # things easier.)
        self.minus_nl2 = MinusInterfaceEnsemble(
            state_vol=vol1,
            innermost_vols=vol1,
            n_l=2
        )
        self.minus_interstitial_nl2 = MinusInterfaceEnsemble(
            state_vol=vol1,
            innermost_vols=vol2,
            n_l=2
        )
        self.minus_nl3 = MinusInterfaceEnsemble(
            state_vol=vol1,
            innermost_vols=vol1,
            n_l=3
        )

    def test_dict_round_trip(self):
        dct = self.minus_nl2.to_dict()
        rebuilt = MinusInterfaceEnsemble.from_dict(dct)
        dct2 = rebuilt.to_dict()
        assert_equal(dct, dct2)

    @raises(ValueError)
    def test_minus_nl1_fail(self):
        minus_nl1 = MinusInterfaceEnsemble(state_vol=vol1,
                                           innermost_vols=vol2,
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
            'lower_in_hit_out_in_out',
            'out_hit_in_out_in'
        ]
        self._test_everything(self.minus_nl2.can_append, non_default, True)

    def test_minus_nl2_strict_can_append(self):
        non_default = [
            'in',
            'lower_hit',
            'in_out',
            'in_out_in',
            'in_out_out_out',
            'in_out_out_in',
            'in_out_in_out',
            'in_out_in_in',
            'in_out_out_in_out',
            'in_out_in_in_out',
            'in_cross',
            'in_cross_in',
            'in_out_cross',
            'in_out_cross_in',
            'in_out_cross_out_in',
            'in_cross_in_cross',
            'upper_in_hit_out_in_out',
            'upper_in_hit_in',
            'upper_in_hit_out'
        ]
        self._test_everything(self.minus_nl2.strict_can_append, non_default,
                              False)

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

    def test_minus_nl2_strict_can_prepend(self):
        non_default = [
            'in',
            'lower_hit',
            'out_in',
            'cross_in',
            'in_cross_in',
            'in_out_cross_in',
            'in_out_cross_out_in',
            'in_out_in',
            'out_out_in',
            'out_out_out_in',
            'in_in_cross_in',
            'out_in_cross_in',
            'out_in_out_out_in',
            'in_in_out_in',
            'in_out_out_in',
            'out_in_out_in',
            'out_in_in_out_in',
            'out_out_in_out_in',
            'cross_in_cross_in',
            'out_hit_in_out_in',
            'upper_in_hit_in',
            'upper_out_hit_in',
        ]
        self._test_everything(self.minus_nl2.strict_can_prepend,
                              non_default, False)

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

    def test_minus_interstitial_nl2_strict_can_append(self):
        non_default = [
            'in',
            'lower_hit',
            'in_out',
            'in_out_out_out',
            'in_cross',
            'in_cross_in',
            'in_out_cross',
            'in_out_cross_in',
            'in_out_cross_out_in',
            'in_cross_in_cross',
            'upper_in_hit_out'
        ]
        self._test_everything(self.minus_interstitial_nl2.strict_can_append,
                              non_default, False)

    def test_minus_interstitial_nl2_can_prepend(self):
        non_default = [
            'in_cross_in_cross_in',
            'in_out_cross_out_in_out_in_out_cross_out_in',
            'in_cross_in_cross'
        ]
        self._test_everything(self.minus_interstitial_nl2.can_prepend,
                              non_default, True)

    def test_minus_interstitial_nl2_strict_can_prepend(self):
        non_default = [
            'in',
            'lower_hit',
            'out_in',
            'cross_in',
            'in_cross_in',
            'in_out_cross_in',
            'in_out_cross_out_in',
            'out_out_in',
            'out_out_out_in',
            'in_in_cross_in',
            'out_in_cross_in',
            'cross_in_cross_in',
            'upper_out_hit_in',
        ]
        self._test_everything(self.minus_interstitial_nl2.strict_can_prepend,
                              non_default, False)

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

    def test_minus_nl3_strict_can_append(self):
        non_default = [
            'in',
            'lower_hit',
            'in_out',
            'in_out_in',
            'in_out_out_out',
            'in_out_out_in',
            'in_out_in_out',
            'in_out_in_in',
            'in_out_out_in_out',
            'in_out_in_in_out',
            'in_cross',
            'in_cross_in',
            'in_out_cross',
            'in_out_cross_in',
            'in_out_cross_out_in',
            'in_cross_in_cross',
            'in_cross_in_cross_in',
            'in_out_in_in_out_in',
            'in_out_in_out_in',
            'upper_in_hit_out_in_out',
            'upper_in_hit_in',
            'upper_in_hit_out',
        ]
        self._test_everything(self.minus_nl3.strict_can_append, non_default,
                              False)

    def test_minus_nl3_can_prepend(self):
        non_default = [
            'in_out_cross_out_in_out_in_out_cross_out_in',
            'out_in_in_in_out_in_out_in_in_in_out'
        ]
        self._test_everything(self.minus_nl3.can_prepend, non_default, True)

    def test_minus_nl3_strict_can_prepend(self):
        non_default = [
            'in',
            'lower_hit',
            'out_in',
            'cross_in',
            'in_cross_in',
            'in_out_cross_in',
            'in_out_cross_out_in',
            'in_out_in',
            'out_out_in',
            'out_out_out_in',
            'in_in_cross_in',
            'out_in_cross_in',
            'out_in_out_out_in',
            'in_in_out_in',
            'in_out_out_in',
            'out_in_out_in',
            'out_in_in_out_in',
            'out_out_in_out_in',
            'cross_in_cross_in',
            'out_hit_in_out_in',
            'in_cross_in_cross_in',
            'in_out_in_in_out_in',
            'in_out_in_out_in',
            'upper_in_hit_in',
            'upper_out_hit_in'
        ]
        self._test_everything(self.minus_nl3.strict_can_prepend,
                              non_default, False)

    def test_extend_sample_from_trajectories(self):
        # set up ensA and ensB
        ensA = paths.TISEnsemble(vol1, vol3, vol1, op)
        ensB = paths.TISEnsemble(vol1, vol3, vol2, op)
        # set up trajA and trajB
        trajA = make_1d_traj([0.25, 1.0, 1.5, 2.1])
        trajB = ttraj['upper_in_cross_in']

        sset = paths.SampleSet([
            paths.Sample(replica=0, ensemble=ensA, trajectory=trajA),
            paths.Sample(replica=1, ensemble=ensB, trajectory=trajB)
        ])
        sset.sanity_check()

        # test with first trajectory
        predestined_snaps = [trajB[-1]]+ttraj['upper_out_in']
        predestined_traj = [s.xyz[0][0] for s in predestined_snaps]
        engine = CalvinistDynamics(predestined_traj)
        sample = self.minus_nl2.extend_sample_from_trajectories(
            sset, replica=-1, engine=engine, level='complex'
        )

        assert_equal(sample.ensemble(sample.trajectory), True)
        assert_equal(sample.ensemble, self.minus_nl2)
        assert_equal(sample.replica, -1)
        assert_equal(len(sample.trajectory), 5)
        expected = trajB + ttraj['upper_out_in']
        for (t, b) in zip(sample.trajectory, expected):
            assert_equal(t.xyz[0][0], b.xyz[0][0])

        # test with a different trajectory
        predestined_snaps = [trajB[-1]]+ttraj['upper_in_out_in']
        predestined_traj = [s.xyz[0][0] for s in predestined_snaps]
        engine = CalvinistDynamics(predestined_traj)
        sample = self.minus_nl2.extend_sample_from_trajectories(
            sset, replica=-1, engine=engine, level='complex'
        )

        assert_equal(sample.ensemble(sample.trajectory), True)
        assert_equal(sample.ensemble, self.minus_nl2)
        assert_equal(sample.replica, -1)
        assert_equal(len(sample.trajectory), 6)
        expected = trajB + ttraj['upper_in_out_in']
        for (t, b) in zip(sample.trajectory, expected):
            assert_equal(t.xyz[0][0], b.xyz[0][0])


# TODO: this whole class should become a single test in SeqEns
class TestSingleEnsembleSequentialEnsemble(EnsembleTest):
    def setup(self):
        #self.inner_ens = AllInXEnsemble(vol1 | vol2)
        self.inner_ens = LengthEnsemble(3) & AllInXEnsemble( vol1 | vol2 )
        self.ens = SequentialEnsemble([self.inner_ens])

    def test_it_all(self):
        for test in list(ttraj.keys()):
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.ens, ttraj[test],
                              self.inner_ens(ttraj[test]), failmsg)
            self._single_test(self.ens.can_append, ttraj[test],
                              self.inner_ens.can_append(ttraj[test]), failmsg)
            self._single_test(self.ens.can_prepend, ttraj[test],
                              self.inner_ens.can_prepend(ttraj[test]), failmsg)


class TestEnsembleSplit(EnsembleTest):
    def setup(self):
        self.inA = AllInXEnsemble(vol1)
        self.outA = AllOutXEnsemble(vol1)

    def test_split(self):
#        raise SkipTest
        traj1 = ttraj['upper_in_out_in_in']
        # print [s for s in traj1]
        subtrajs_in_1 = self.inA.split(traj1)
        # print subtrajs_in_1
        # print [[s for s in t] for t in subtrajs_in_1]
        assert_equal(len(subtrajs_in_1), 2)
        assert_equal(len(subtrajs_in_1[0]), 1)
        assert_equal(len(subtrajs_in_1[1]), 2)
        subtrajs_out_1 = self.outA.split(traj1)
        assert_equal(len(subtrajs_out_1), 1)

        traj2 = ttraj['upper_in_out_in_in_out_in']
        # print [s for s in traj2]
        subtrajs_in_2 = self.inA.split(traj2)
        # print [[s for s in t] for t in subtrajs_in_2]
        assert_equal(len(subtrajs_in_2), 3)
        assert_equal(len(subtrajs_in_2[0]), 1)
        assert_equal(len(subtrajs_in_2[1]), 2)
        assert_equal(len(subtrajs_in_2[2]), 1)
        subtrajs_out_2 = self.outA.split(traj2)
        assert_equal(len(subtrajs_out_2), 2)
        assert_equal(len(subtrajs_out_2[0]), 1)
        assert_equal(len(subtrajs_out_2[1]), 1)

        ensembleAXA = paths.SequentialEnsemble([
            self.inA,
            self.outA,
            self.inA
        ])

        traj3 = make_1d_traj(coordinates=[0.3, 0.6, 0.3, 0.6, 0.3])

        assert(self.inA(paths.Trajectory([traj3[0]])))
        assert(self.outA(paths.Trajectory([traj3[1]])))

        subtrajs_in_3 = ensembleAXA.split(traj3)
        assert_equal((len(subtrajs_in_3)), 2)
        assert_equal((len(subtrajs_in_3[0])), 3)
        assert_equal((len(subtrajs_in_3[1])), 3)
        assert(traj3.subtrajectory_indices(subtrajs_in_3[0]) == [0, 1, 2])
        assert(traj3.subtrajectory_indices(subtrajs_in_3[1]) == [2, 3, 4])

        subtrajs_in_3 = ensembleAXA.split(traj3, reverse=True)
        assert_equal((len(subtrajs_in_3)), 2)
        assert_equal((len(subtrajs_in_3[0])), 3)
        assert_equal((len(subtrajs_in_3[1])), 3)
        assert(traj3.subtrajectory_indices(subtrajs_in_3[0]) == [2, 3, 4])
        assert(traj3.subtrajectory_indices(subtrajs_in_3[1]) == [0, 1, 2])

        subtrajs_in_3 = ensembleAXA.split(traj3, overlap=0)
        assert_equal((len(subtrajs_in_3)), 1)
        assert_equal((len(subtrajs_in_3[0])), 3)
        assert(traj3.subtrajectory_indices(subtrajs_in_3[0]) == [0, 1, 2])

        subtrajs_in_3 = ensembleAXA.split(traj3, reverse=True, overlap=0)
        assert_equal((len(subtrajs_in_3)), 1)
        assert_equal((len(subtrajs_in_3[0])), 3)
        assert(traj3.subtrajectory_indices(subtrajs_in_3[0]) == [2, 3, 4])

        subtrajs_in_3 = ensembleAXA.split(traj3, overlap=1, max_length=2)
        assert_equal((len(subtrajs_in_3)), 0)

        subtrajs_in_3 = ensembleAXA.split(traj3, reverse=True, max_length=2)
        assert_equal((len(subtrajs_in_3)), 0)

        subtrajs_in_3 = ensembleAXA.split(traj3, max_length=3)
        assert_equal(len(subtrajs_in_3), 2)
        assert_equal((len(subtrajs_in_3[0])), 3)
        assert_equal((len(subtrajs_in_3[1])), 3)
        assert(traj3.subtrajectory_indices(subtrajs_in_3[0]) == [0, 1, 2])
        assert(traj3.subtrajectory_indices(subtrajs_in_3[1]) == [2, 3, 4])

        subtrajs_in_3 = ensembleAXA.split(traj3, reverse=True, max_length=3)
        assert_equal((len(subtrajs_in_3)), 2)
        assert_equal((len(subtrajs_in_3[0])), 3)
        assert_equal((len(subtrajs_in_3[1])), 3)
        assert(traj3.subtrajectory_indices(subtrajs_in_3[1]) == [0, 1, 2])
        assert(traj3.subtrajectory_indices(subtrajs_in_3[0]) == [2, 3, 4])

        subtrajs_in_3 = ensembleAXA.split(traj3, reverse=False, min_length=4)
        assert_equal((len(subtrajs_in_3)), 0)

        subtrajs_in_3 = ensembleAXA.split(traj3, reverse=True, min_length=4)
        assert_equal((len(subtrajs_in_3)), 0)

        sub_traj = ensembleAXA.find_first_subtrajectory(traj3)
        assert(traj3.subtrajectory_indices(sub_traj) == [0,1,2])

        sub_traj = ensembleAXA.find_last_subtrajectory(traj3)
        assert(traj3.subtrajectory_indices(sub_traj) == [2,3,4])

class TestVolumeCombinations(EnsembleTest):
    def setup(self):
        self.outA = paths.AllOutXEnsemble(vol1)
        self.outB = paths.AllOutXEnsemble(~vol2)
        self.outA.special_debug = True
        self.outB.special_debug = True
        self.partinA = paths.PartInXEnsemble(vol1)
        self.partinB = paths.PartInXEnsemble(~vol2)
        self.outA_or_outB = self.outA | self.outB
        self.outA_and_outB = self.outA & self.outB
        self.partinA_or_partinB = self.partinA | self.partinB
        self.partinA_and_partinB = self.partinA & self.partinB
        extras = build_trajdict(['babbc', 'ca', 'bcbba', 'abbc', 'cbba',
                                 'abbcb', 'cbbab'], lower, upper)
        for test in list(extras.keys()):
            extras[test] = make_1d_traj(coordinates=extras[test],
                                       velocities=[1.0]*len(extras[test]))
        self.local_ttraj = dict(ttraj)
        self.local_ttraj.update(extras)

    def _test_trusted(self, trajectory, function, results,
                      cache_results=None, direction=+1, start_traj_len=1):
        # Tests `trajectory` frame by frame in a forward direction for the
        # `function`, expecting `results`. Additionally, can take the 

        if cache_results is None:
            cache_results = {}

        # clear the caches before starting
        for cache in list(cache_results.keys()):
            cache.__init__(direction=cache.direction)

        for i in range(len(trajectory)-start_traj_len):
            if direction > 0:
                start = 0
                end = start + (i+start_traj_len)
            elif direction < 0:
                end = len(trajectory)
                start = end - (i+start_traj_len)
            # test untrusted
            assert_equal(function(trajectory[start:end]), results[i])
            # test trusted
            trusted_val = function(trajectory[start:end], trusted=True)
            # print i, "["+str(start)+":"+str(end)+"]", trusted_val, results[i]
            assert_equal(trusted_val, results[i])
            for cache in list(cache_results.keys()):
                # TODO: this is currently very specific to the caches used
                # by volumes ensembles. That should be generalized by
                # allowing several different tags within contents.
                # cache_results could {cache : {'content_key' : [values]}}
                if cache_results[cache][i] is not None:
                    #print "cache", cache_results.keys().index(cache),
                    try:
                        contents = cache.contents['previous']
                    except KeyError:
                        contents = None
                    #print contents, cache_results[cache][i]

                    assert_equal(cache.contents['previous'],
                                 cache_results[cache][i])

    def test_call_outA_or_outB(self):
        # print self.local_ttraj['upper_out_in_out_out_cross'].xyz[:,0,0]
        self._test_trusted(
            trajectory=self.local_ttraj['upper_out_in_out_out_cross'],
            function=self.outA_or_outB,
            results=[True, True, True, True, False],
            cache_results={
                self.outA._cache_call : [True, False, False, False, False],
                self.outB._cache_call : [None, True, True, True, False]
            }
        )
        self._test_trusted(
            trajectory=self.local_ttraj['upper_out_cross_out_out_in'],
            function=self.outA_or_outB,
            results=[True, True, True, True, False],
            cache_results={
                self.outA._cache_call : [True, True, True, True, False],
                self.outB._cache_call : [None, None, None, None, False]
            }
        )
        self._test_trusted(
            trajectory=self.local_ttraj['upper_in_cross'],
            function=self.outA_or_outB,
            results=[True, False],
            cache_results={
                self.outA._cache_call : [False, False],
                self.outB._cache_call : [True, False]
            }
        )
        self._test_trusted(
            trajectory=self.local_ttraj['upper_cross_in'],
            function=self.outA_or_outB,
            results=[True, False],
            cache_results={
                self.outA._cache_call : [True, False],
                self.outB._cache_call : [None, False]
            }
        )

    def test_call_outA_and_outB(self):
        self._test_trusted(
            trajectory=self.local_ttraj['upper_out_in_out_out_cross'],
            function=self.outA_and_outB,
            results=[True, False, False, False, False],
            cache_results={
                # cache for A gets checked first: value of cache for B
                # doesn't matter once cache for A is False (short-circuit)
                self.outA._cache_call : [True, False, False, False, False],
                self.outB._cache_call : [True, None, None, None, None]
            }
        )
        self._test_trusted(
            trajectory=self.local_ttraj['upper_out_cross_out_out_in'],
            function=self.outA_and_outB,
            results=[True, False, False, False, False],
            cache_results={
                self.outA._cache_call : [True, True, True, True, False],
                self.outB._cache_call : [True, False, False, False, None]
            }
        )

    def test_can_append_outA_or_outB(self):
        self._test_trusted(
            trajectory=self.local_ttraj['upper_out_in_out_out_cross'],
            function=self.outA_or_outB.can_append, 
            results=[True, True, True, True, False], 
            cache_results={
                self.outA._cache_can_append : [True, False, False, False, False],
                self.outB._cache_can_append : [None, True, True, True, False]
            }
        )
        self._test_trusted(
            trajectory=self.local_ttraj['upper_out_cross_out_out_in'],
            function=self.outA_or_outB.can_append,
            results=[True, True, True, True, False],
            cache_results={
                self.outA._cache_can_append : [True, True, True, True, False],
                self.outB._cache_can_append : [None, None, None, None, False]
            }
        )
        self._test_trusted(
            trajectory=self.local_ttraj['upper_in_cross'],
            function=self.outA_or_outB.can_append,
            results=[True, False],
            cache_results={
                self.outA._cache_can_append : [False, False],
                self.outB._cache_can_append : [True, False]
            }
        )
        self._test_trusted(
            trajectory=self.local_ttraj['upper_cross_in'],
            function=self.outA_or_outB.can_append,
            results=[True, False],
            cache_results={
                self.outA._cache_can_append : [True, False],
                self.outB._cache_can_append : [None, False]
            }
        )

    def test_call_start_from_later_frame(self):
        self._test_trusted(
            trajectory=self.local_ttraj['upper_out_cross_out_out_in'],
            function=self.outA_or_outB.can_append,
            results=[True, True, True, False],
            cache_results={
                self.outA._cache_can_append : [True, True, True, False],
                self.outB._cache_can_append : [None, None, None, False]
            },
            start_traj_len=2
        )

    def test_can_append_outA_and_outB(self):
        self._test_trusted(
            trajectory=self.local_ttraj['upper_out_in_out_out_cross'],
            function=self.outA_and_outB.can_append,
            results=[True, False, False, False, False],
            cache_results={
                # cache for A gets checked first: value of cache for B
                # doesn't matter once cache for A is False (short-circuit)
                self.outA._cache_can_append : [True, False, False, False, False],
                self.outB._cache_can_append : [True, None, None, None, None]
            }
        )
        self._test_trusted(
            trajectory=self.local_ttraj['upper_out_cross_out_out_in'],
            function=self.outA_and_outB.can_append,
            results=[True, False, False, False, False],
            cache_results={
                self.outA._cache_can_append : [True, True, True, True, False],
                self.outB._cache_can_append : [True, False, False, False, None]
            }
        )

    def test_can_prepend_outA_or_outB(self):
        self._test_trusted(
            trajectory=self.local_ttraj['upper_in_out_out_cross'],
            function=self.outA_or_outB.can_prepend,
            results=[True, True, True, False],
            cache_results={
                self.outA._cache_can_prepend : [True, True, True, False],
                self.outB._cache_can_prepend : [None, None, None, False]
            },
            direction=-1
        )
        self._test_trusted(
            trajectory=self.local_ttraj['upper_cross_out_out_in'],
            function=self.outA_or_outB.can_prepend,
            results=[True, True, True, False],
            cache_results={
                self.outA._cache_can_prepend : [False, False, False, False],
                self.outB._cache_can_prepend : [True, True, True, False]
            },
            direction=-1
        )

    def test_can_prepend_start_from_later_frame(self):
        self._test_trusted(
            trajectory=self.local_ttraj['upper_in_out_out_cross'],
            function=self.outA_or_outB.can_prepend,
            results=[True, True, False],
            cache_results={
                self.outA._cache_can_prepend : [True, True, False],
                self.outB._cache_can_prepend : [None, None, False]
            },
            direction=-1,
            start_traj_len=2
        )


    def test_can_prepend_outA_and_outB(self):
        self._test_trusted(
            trajectory=self.local_ttraj['upper_in_out_out_cross_out'],
            function=self.outA_and_outB.can_prepend,
            results=[True, False, False, False, False],
            cache_results={
                self.outA._cache_can_prepend : [True, True, True, True, False],
                self.outB._cache_can_prepend : [True, False, False, False, False]
            },
            direction=-1
        )
        self._test_trusted(
            trajectory=self.local_ttraj['upper_cross_out_out_in_out'],
            function=self.outA_and_outB.can_prepend,
            results=[True, False, False, False, False],
            cache_results={
                self.outA._cache_can_prepend : [True, False, False, False, False],
                self.outB._cache_can_prepend : [True, None, None, None, None]
            },
            direction=-1
        )


class TestAbstract(object):
    @raises_with_message_like(TypeError, "Can't instantiate abstract class")
    def test_abstract_ensemble(self):
        mover = paths.Ensemble()

    @raises_with_message_like(TypeError, "Can't instantiate abstract class")
    def test_abstract_volumeensemble(self):
        mover = paths.VolumeEnsemble()

class TestEnsembleEquality(object):
    # generic tests for ensemble equality; we use the EmptyEnsemble as
    # example. See:
    # * https://github.com/openpathsampling/openpathsampling/issues/700
    # * https://github.com/openpathsampling/openpathsampling/issues/701
    def test_empty_ensemble_equality(self):
        ens1 = paths.EmptyEnsemble()
        ens2 = paths.EmptyEnsemble()
        assert_true(ens1 == ens2)
        assert_false(ens1 != ens2)

    # TODO: may add tests for other ensembles, or may move this test
    # somewhere else
