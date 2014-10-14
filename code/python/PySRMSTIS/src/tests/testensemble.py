from nose.tools import assert_equal, assert_not_equal, raises
from duckpunching import CallIdentity, prepend_exception_message

import os
import sys
sys.path.append(os.path.abspath('../'))
from volume import LambdaVolume
from ensemble import *

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
    ''' Setup for tests tests of classes in ensemble.py. '''
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

    crosses = [-0.05, 0.05, 0.1]
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

class EnsembleTest(object):
    def _single_test(self, ensemble, traj, res, failmsg):
        try:
            assert_equal(res, ensemble(traj))
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
        # longest ttraj is 5 = 8-3 frames long
        self.slice_ens = ExitsXEnsemble(vol1, slice(3,8))
        self.wrapstart = 3
        self.wrapend = 10

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
        # longest ttraj is 5 = 8-3 frames long
        self.slice_ens = EntersXEnsemble(vol1, slice(3,8))
        self.wrapstart = 3
        self.wrapend = 10

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
        self.inX = InXEnsemble(vol1)
        self.outX = OutXEnsemble(vol1)
        self.hitX = HitXEnsemble(vol1)
        self.leaveX = LeaveXEnsemble(vol1)
        self.enterX = EntersXEnsemble(vol1)
        self.exitX = ExitsXEnsemble(vol1)

    @raises(ValueError)
    def test_maxminoverlap_size(self):
        """SequentialEnsemble error if max/min overlap sizes different"""
        SequentialEnsemble([self.inX, self.outX, self.inX], (0,0), (0,0,0))

    @raises(ValueError)
    def test_maxoverlap_ensemble_size(self):
        """SequentialEnsemble error if overlap sizes don't match ensemble size"""
        SequentialEnsemble([self.inX, self.outX, self.inX], (0,0,0), (0,0,0))

    @raises(ValueError)
    def test_minmax_order(self):
        """SequentialEnsemble error if min_overlap > max_overlap"""
        SequentialEnsemble([self.inX, self.outX, self.inX], (0,1), (0,0))

    def test_allowed_initializations(self):
        """SequentialEnsemble initializes correctly"""
        A = SequentialEnsemble([self.inX, self.outX, self.inX], (0,0), (0,0))
        B = SequentialEnsemble([self.inX, self.outX, self.inX],0,0)
        C = SequentialEnsemble([self.inX, self.outX, self.inX])
        assert_equal(A.min_overlap,B.min_overlap)
        assert_equal(A.min_overlap,C.min_overlap)
        assert_equal(A.max_overlap,B.max_overlap)
        assert_equal(A.max_overlap,C.max_overlap)

    def test_sequential_in_out(self):
        """SequentialEnsembles based on In/OutXEnsemble"""
        # idea: for each ttraj, use the key name to define in/out behavior,
        # dynamically construct a 
        pass

    def test_sequential_hit_leave(self):
        """SequentialEnsembles based on Hit/LeaveXEnsemble"""
        pass

    def test_sequential_enter_exit(self):
        """SequentialEnsembles based on Enters/ExitsXEnsemble"""
        pass

    def test_sequential_str(self):
        pass
