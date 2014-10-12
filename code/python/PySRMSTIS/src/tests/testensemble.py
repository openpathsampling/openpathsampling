from nose.tools import assert_equal, assert_not_equal, raises
from duckpunching import CallIdentity

import os
import sys
sys.path.append(os.path.abspath('../'))
from volume import LambdaVolume
from ensemble import ExitsXEnsemble, EntersXEnsemble

def wrap_traj(traj, slice0, sliceN):
    """Wraps the traj such that the original traj is given by
    traj[slice(slice0,sliceN)], by padding beginning with traj[0] and
    end with traj[-1]. Used to test the slice restricted
    trajectories."""
    #TODO
    pass

def setUp():
    ''' Setup for tests tests of classes in ensemble.py. '''
    lower = 0.1
    upper = 0.5
    global op = CallIdentity()

    global vol1 = LambdaVolume(op, lower, upper)

    all_positive = [0.1, 0.05, 0.15]
    # generate trajectories always outside the border
    global traj_lower_in = map(lower.__add__, all_positive)
    global traj_upper_in = map(upper.__sub__, all_positive)
    global traj_lower_out = map(lower.__sub__, all_positive)
    global traj_upper_out = map(upper.__add__, all_positive)

    hits_border = [0.1, 0.0, 0.15]
    # generate trajectories which hit the border but do not cross
    global traj_lower_in_hit_in = map(lower.__add__, hits_border)
    global traj_lower_out_hit_out = map(lower.__sub__, hits_border)
    global traj_upper_in_hit_in = map(upper.__sub__, hits_border)
    global traj_upper_out_hit_out = map(upper.__add__, hits_border)

    crosses = [-0.05, 0.05, 0.1]
    border_then_cross = [-0.1, 0.0, 0.1]
    # generate trajectories which cross the border from the inside
    global traj_lower_in_out = map(lower.__sub__, crosses)
    global traj_lower_in_hit_out = map(lower.__sub__, border_then_cross)
    global traj_upper_in_out = map(upper.__add__, crosses)
    global traj_upper_in_hit_out = map(upper.__add__, border_then_cross)

    # generate trajectories which cross the border from the outside
    global traj_lower_out_in = map(lower.__add__, crosses)
    global traj_lower_out_hit_in = map(lower.__add__, border_then_cross)
    global traj_upper_out_in = map(upper.__sub__, crosses)
    global traj_upper_out_hit_in = map(upper.__sub__, border_then_cross)

    # generate trajectories which cross the border in both directions
    #TODO


class testExitsXEnsemble(object):
    def setUp(self):
        self.fulltraj_ensemble = ExitsXEnsemble(vol1)
        self.parttraj_ensemble = ExitsXEnsemble(vol1, slice(2,6))

    @raises(Exception)
    def test_single_frame(self):
        '''ExitsXEnsemble built with single frame raises Exception'''
        single_frame_ensemble = ExitsXEnsemble(vol1, 2)

    @raises(Exception)
    def test_modified_to_single_frame(self)
        '''ExitsXEnsemble modified to have single frame raises Exception'''
        single_frame_ensemble = ExitsXEnsemble(vol1)
        single_frame_ensemble.frames = 2
        single_frame_ensemble(traj_in_out_upper)

    def test_noncrossing(self):
        pass

    def test_exit(self):
        pass

    def test_entrance(self):
        pass

    def test_str(self):
        pass


class testEntersXEnsemble(object):
    def setUp(self):
        pass

    def teardown(self):
        pass

    def test_single_frame(self):
        pass

    def test_noncrossing(self):
        pass

    def test_exit(self):
        pass

    def test_entrance(self):
        pass

    def test_str(self):
        pass
