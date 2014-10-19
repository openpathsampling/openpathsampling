from nose.tools import assert_equal, assert_not_equal, raises
from nose.plugins.skip import Skip, SkipTest
from duckpunching import CallIdentity, prepend_exception_message

import os
import sys
sys.path.append(os.path.abspath('../'))
from trajectory import Trajectory
from volume import LambdaVolume
from ensemble import *

import re

def wrap_traj(traj, start, length):
    """Wraps the traj such that the original traj starts at frame `start`
    and is of length `length` by padding beginning with traj[0] and end with
    traj[-1]. Used to test the slice restricted trajectories."""
    if (start < 0) or (length < len(traj)+start):
        raise ValueError("""wrap_traj: start < 0 or length < len(traj)+start
                                {0} < 0 or {1} < {2}+{0}""".format(
                                start, length, len(traj)) )
    outtraj = traj[:] # shallow copy
    # prepend
    for i in range(start):
        outtraj.insert(0,traj[0])
    # append
    for i in range(length - (len(traj)+start)):
        outtraj.append(traj[-1])
    return outtraj

def test_wrap_traj():
    """Testing wrap_traj (oh gods, the meta! a test for a test function!)"""
    intraj = [1,2,3]
    assert_equal(wrap_traj(intraj, 3, 6), [1, 1, 1, 1, 2, 3])
    assert_equal(wrap_traj(intraj, 3, 8), [1, 1, 1, 1, 2, 3, 3, 3])
    assert_equal(wrap_traj(intraj, 3, 8)[slice(3, 6)], intraj)

def setUp():
    ''' Setup for tests of classes in ensemble.py. '''
    global lower, upper, op, vol1, ttraj
    lower = 0.1
    upper = 0.5
    op = CallIdentity()
    vol1 = LambdaVolume(op, lower, upper)
    ttraj = {}

    all_positive = [0.1, 0.05, 0.15]
    # generate trajectories always outside the border
    ttraj['lower_in'] = map(lower.__add__, all_positive)
    ttraj['upper_in'] = map(upper.__sub__, all_positive)
    ttraj['lower_out'] = map(lower.__sub__, all_positive)
    ttraj['upper_out'] = map(upper.__add__, all_positive)

    hits_border = [0.1, 0.0, 0.15]
    # generate trajectories which hit the border but do not cross
    ttraj['lower_in_hit_in'] = map(lower.__add__, hits_border)
    ttraj['lower_out_hit_out'] = map(lower.__sub__, hits_border)
    ttraj['upper_in_hit_in'] = map(upper.__sub__, hits_border)
    ttraj['upper_out_hit_out'] = map(upper.__add__, hits_border)

    crosses = [-0.05, 0.05]
    border_then_cross = [-0.1, 0.0, 0.1]
    # generate trajectories which cross the border from the inside
    ttraj['lower_in_out'] = map(lower.__sub__, crosses)
    ttraj['lower_in_hit_out'] = map(lower.__sub__, border_then_cross)
    ttraj['upper_in_out'] = map(upper.__add__, crosses)
    ttraj['upper_in_hit_out'] = map(upper.__add__, border_then_cross)

    # generate trajectories which cross the border from the outside
    ttraj['lower_out_in'] = map(lower.__add__, crosses)
    ttraj['lower_out_hit_in'] = map(lower.__add__, border_then_cross)
    ttraj['upper_out_in'] = map(upper.__sub__, crosses)
    ttraj['upper_out_hit_in'] = map(upper.__sub__, border_then_cross)

    doublecross = [-0.1, 0.0, 0.1, -0.05, 0.05]
    # generate trajectories which cross the border in both directions
    ttraj['lower_out_in_out_in'] = map(lower.__add__, doublecross)
    ttraj['lower_in_out_in_out'] = map(lower.__sub__, doublecross)
    ttraj['upper_out_in_out_in'] = map(upper.__sub__, doublecross)
    ttraj['upper_in_out_in_out'] = map(upper.__add__, doublecross)

    ttraj['lower_hit_1'] = [lower]
    ttraj['upper_hit_1'] = [upper]

    aba = [-0.05, 0.05, -0.1]
    ttraj['lower_in_out_in'] = map(lower.__sub__, aba)
    ttraj['upper_in_out_in'] = map(upper.__add__, aba)
    ttraj['lower_out_in_out'] = map(lower.__add__, aba)
    ttraj['upper_out_in_out'] = map(upper.__sub__, aba)
    abaab = [-0.05, 0.05, -0.1, -0.05, 0.1]
    a = [-0.05]
    ttraj['upper_in_1'] = map(upper.__add__, a)
    ttraj['lower_in_1'] = map(lower.__sub__, a)
    ttraj['upper_out_1'] = map(upper.__sub__, a)
    ttraj['lower_out_1'] = map(lower.__add__, a)
    abaa = [-0.05, 0.05, -0.1, -0.05]
    ttraj['upper_in_out_in_in'] = map(upper.__add__, abaa)
    ttraj['lower_in_out_in_in'] = map(lower.__sub__, abaa)
    ttraj['upper_out_in_out_out'] = map(upper.__sub__, abaa)
    ttraj['lower_out_in_out_out'] = map(lower.__add__, abaa)
    ababa = [-0.05, 0.05, -0.1, 0.1, -0.05]
    ttraj['upper_in_out_in_out_in'] = map(upper.__add__, ababa)
    ttraj['lower_in_out_in_out_in'] = map(lower.__sub__, ababa)
    ttraj['upper_out_in_out_in_out'] = map(upper.__sub__, ababa)
    ttraj['lower_out_in_out_in_out'] = map(lower.__add__, ababa)
    abaaba = [-0.05, 0.05, -0.1, -0.15, 0.05, -0.05]
    ttraj['upper_in_out_in_in_out_in'] = map(upper.__add__, abaaba)
    ttraj['lower_in_out_in_in_out_in'] = map(lower.__sub__, abaaba)
    ttraj['upper_out_in_out_out_in_out'] = map(upper.__sub__, abaaba)
    ttraj['lower_out_in_out_out_in_out'] = map(lower.__add__, abaaba)
    abaab = [-0.05, 0.05, -0.15, -0.1, 0.05]
    ttraj['upper_in_out_in_in_out'] = map(upper.__add__, abaab)
    ttraj['lower_in_out_in_in_out'] = map(lower.__sub__, abaab)
    ttraj['upper_out_in_out_out_in'] = map(upper.__sub__, abaab)
    ttraj['lower_out_in_out_out_in'] = map(lower.__add__, abaab)
    abba = [-0.05, 0.1, 0.05, -0.1]
    ttraj['upper_in_out_out_in'] = map(upper.__add__, abba)
    ttraj['lower_in_out_out_in'] = map(lower.__sub__, abba)
    ttraj['upper_out_in_in_out'] = map(upper.__sub__, abba)
    ttraj['lower_out_in_in_out'] = map(lower.__add__, abba)
    abbab = [-0.05, 0.1, 0.05, -0.1, 0.15]
    ttraj['upper_in_out_out_in_out'] = map(upper.__add__, abbab)
    ttraj['lower_in_out_out_in_out'] = map(lower.__sub__, abbab)
    ttraj['upper_out_in_in_out_in'] = map(upper.__sub__, abbab)
    ttraj['lower_out_in_in_out_in'] = map(lower.__add__, abbab)

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
        if to_append != None:
            if res == []:
                res.append(to_append)
            elif to_append != res[-1]:
                res.append(to_append)
    return res

class EnsembleTest(object):
    def _single_test(self, ensemble_fcn, traj, res, failmsg):
        try:
            assert_equal(res, ensemble_fcn(traj))
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
        results = { 'upper_in' : False,
                    'upper_out' : False,
                    'lower_in' : False,
                    'lower_out' : False
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
        self.length1 = LengthEnsemble(1)
        # properly speaking, tis and minus ensembles require outX and
        # leave_interface; but for now I'll just do the easy way
        self.tis_ensemble = SequentialEnsemble( [
                                    self.inX & self.length1,
                                    self.outX,
                                    self.inX & self.length1 ]
                                    )
        self.minus_ensemble = SequentialEnsemble( [
                                    self.inX & self.length1,
                                    self.outX,
                                    self.inX,
                                    self.outX,
                                    self.inX & self.length1 ]
                                    )

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
                        'upper_in_1' : True,
                        'lower_in_1' : True,
                        'upper_in' : False,
                        'lower_in' : False,
                        'upper_out' : True,
                        'lower_out' : True,
                        'upper_out_in' : False,
                        'lower_out_in' : False,
                        'upper_out_1' : True,
                        'lower_out_1' : True,
                        'upper_in_out_in_in' : False,
                        'lower_in_out_in_in' : False,
                        'upper_in_out_in_out_in' : False,
                        'lower_in_out_in_out_in' : False
                    }   
        for test in results.keys():
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.tis_ensemble.forward, 
                                ttraj[test], results[test], failmsg)

    def test_can_append_minus(self):
        """SequentialEnsemble as MinusEnsemble knows when it can append"""
        results =   {   'upper_in_out' : True,
                        'lower_in_out' : True,
                        'upper_in_out_in' : True,
                        'lower_in_out_in' : True,
                        'upper_in_1' : True,
                        'lower_in_1' : True,
                        'upper_in' : True,
                        'lower_in' : True,
                        'upper_out' : True,
                        'lower_out' : True,
                        'upper_out_in' : True,
                        'lower_out_in' : True,
                        'upper_out_1' : True,
                        'lower_out_1' : True,

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
            self._single_test(self.minus_ensemble.forward, 
                                ttraj[test], results[test], failmsg)

    def test_can_prepend_tis(self):
        """SequentialEnsemble as TISEnsemble knows when it can prepend"""
        results =   {   'upper_in_out' : False,
                        'lower_in_out' : False,
                        'upper_in_out_in' : False,
                        'lower_in_out_in' : False,
                        'upper_in_1' : True,
                        'lower_in_1' : True,
                        'upper_in' : False,
                        'lower_in' : False,
                        'upper_out' : True,
                        'lower_out' : True,
                        'upper_out_in' : True,
                        'lower_out_in' : True,
                        'upper_out_1' : True,
                        'lower_out_1' : True,
                        'upper_in_out_in_in' : False,
                        'lower_in_out_in_in' : False,
                        'upper_in_out_in_out_in' : False,
                        'lower_in_out_in_out_in' : False
                    }   
        for test in results.keys():
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.tis_ensemble.backward, 
                                ttraj[test], results[test], failmsg)


    def test_can_prepend_minus(self):
        """SequentialEnsemble as MinusEnsemble knows when it can prepend"""
        results =   {   'upper_in_out' : True,
                        'lower_in_out' : True,
                        'upper_in_out_in' : True,
                        'lower_in_out_in' : True,
                        'upper_in_1' : True,
                        'lower_in_1' : True,
                        'upper_in' : True,
                        'lower_in' : True,
                        'upper_out' : True,
                        'lower_out' : True,
                        'upper_out_in' : True,
                        'lower_out_in' : True,
                        'upper_out_1' : True,
                        'lower_out_1' : True,

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
            self._single_test(self.minus_ensemble.backward, 
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

    def test_sequential_tis(self):
        """SequentialEnsemble as TISEnsemble identifies paths"""
        results = {}
        for test in ttraj.keys():
            results[test] = False
        results['upper_in_out_in'] = True
        results['lower_in_out_in'] = True
        results['upper_in_out_out_in'] = True
        results['lower_in_out_out_in'] = True
        for test in results.keys():
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.tis_ensemble, ttraj[test], results[test],
                              failmsg)

    def test_sequential_minus(self):
        """SequentialEnsemble as MinusEnsemble identifies paths"""
        results = {}
        for test in ttraj.keys():
            results[test] = False
        results['upper_in_out_in_out_in'] = True
        results['lower_in_out_in_out_in'] = True
        results['upper_in_out_in_in_out_in'] = True
        results['lower_in_out_in_in_out_in'] = True
        for test in results.keys():
            failmsg = "Failure in "+test+"("+str(ttraj[test])+"): "
            self._single_test(self.minus_ensemble, ttraj[test], results[test],
                              failmsg)

    def test_sequential_enter_exit(self):
        """SequentialEnsembles based on Enters/ExitsXEnsemble"""
        # TODO: this includes a test of the overlap ability
        raise SkipTest

    def test_str(self):
        #TODO: passes test, but isn't right -- need to fix another part of
        #string function somewhere
        assert_equal(self.tis_ensemble.__str__(), """[
(
  x[t] in {x|Id(x) in [0.1, 0.5]} for all t in [None:None])
)
and
(
  len(x) = 1
),
x[t] not in (not {x|Id(x) in [0.1, 0.5]}) for all t in [None:None]),
(
  x[t] in {x|Id(x) in [0.1, 0.5]} for all t in [None:None])
)
and
(
  len(x) = 1
)
]""")
        self.minus_ensemble.__str__()
