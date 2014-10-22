"""
@author David W.H. Swenson
"""

import os
import sys
from nose.tools import assert_equal, assert_not_equal, raises
from nose.plugins.skip import Skip, SkipTest
from test_helpers import CallIdentity

sys.path.append(os.path.abspath('../'))
import volume

def setUp():
    global op_id, lambdavol
    op_id = CallIdentity()
    lambdavol = volume.LambdaVolume(op_id, -0.5, 0.5)


class testEmptyVolume(object):
    def test_empty_volume(self):
        """Empty volume is well-behaved"""
        empty = volume.EmptyVolume()
        test = 0.1
        assert_equal(empty(test), False)
        assert_equal((empty & lambdavol)(test), False)
        assert_equal((lambdavol & empty)(test), False)
        assert_equal((empty | lambdavol)(test), True)
        assert_equal((lambdavol | empty)(test), True)
        assert_equal((empty & lambdavol).__str__(), "empty")
        assert_equal((~ empty).__str__(), "all")

class testFullVolume(object):
    def test_full_volume(self):
        """Full volume is well-behaved"""
        full = volume.FullVolume()
        test = 0.1
        assert_equal(full(test), True)
        assert_equal((full & lambdavol)(test), True)
        assert_equal((lambdavol & full)(test), True)
        assert_equal((full | lambdavol)(test), True)
        assert_equal((lambdavol | full)(test), True)
        assert_equal((lambdavol | full).__str__(), "all")

class testLambdaVolume(object):
    pass

class testLambdaVolumePeriodic(object):
    
    def test_normal(self):
        """min<max and both within periodic domain"""
        lambda_min = -150
        lambda_max = 70
        vol = volume.LambdaVolumePeriodic(op_id,
                                          lambda_min, lambda_max, -180,180)
        assert_equal(vol.period_len, 360)
        assert_equal(vol.period_shift, -180)
        assert_equal(vol.__str__(), 
            "{x|(Id(x) - -180) % 360 + -180 in [-150, 70]}")
        for periodic_image in [-1, 0, 1]:
            # out of state
            assert_equal(False, vol(lambda_min-1.0 + 360*periodic_image))
            assert_equal(False, vol(lambda_max+1.0 + 360*periodic_image))
            # in state
            assert_equal(True, vol(lambda_min+1.0 + 360*periodic_image))
            assert_equal(True, vol(lambda_max-1.0 + 360*periodic_image))
            # border
            # TODO: consider changing upper border from True to False
            assert_equal(True, vol(lambda_min + 360*periodic_image))
            assert_equal(True, vol(lambda_max + 360*periodic_image))

    def test_normal_implicit(self):
        """min<max, no periodic domain defined"""
        lambda_min = -150
        lambda_max = 70
        vol = volume.LambdaVolumePeriodic(op_id,
                                          lambda_min, lambda_max)
        assert_equal(vol.__str__(),
            "{x|Id(x) [periodic] in [-150, 70]}")
        # out of state
        assert_equal(False, vol(lambda_min-1.0))
        assert_equal(False, vol(lambda_max+1.0))
        # in state
        assert_equal(True, vol(lambda_min+1.0))
        assert_equal(True, vol(lambda_max-1.0))
        # border
        assert_equal(True, vol(lambda_min))
        assert_equal(True, vol(lambda_max))

    def test_inverted(self):
        """max<min and both within periodic domain"""
        lambda_min = 70
        lambda_max = -150
        vol = volume.LambdaVolumePeriodic(op_id,
                                          lambda_min, lambda_max, -180,180)
        assert_equal(vol.__str__(),
            "{x|(Id(x) - -180) % 360 + -180 in [-180, -150] union [70, 180]}")
        for periodic_image in [-1, 0, 1]:
            # out of state
            assert_equal(False, vol(lambda_min-1.0 + 360*periodic_image))
            assert_equal(False, vol(lambda_max+1.0 + 360*periodic_image))
            # in state
            assert_equal(True, vol(lambda_min+1.0 + 360*periodic_image))
            assert_equal(True, vol(lambda_max-1.0 + 360*periodic_image))
            # border
            assert_equal(True, vol(lambda_min + 360*periodic_image))
            assert_equal(True, vol(lambda_max + 360*periodic_image))

    def test_inverted_implicit(self):
        """max<min, no periodic domain defined"""
        lambda_min = 70
        lambda_max = -150
        vol = volume.LambdaVolumePeriodic(op_id,
                                          lambda_min, lambda_max)
        # out of state
        assert_equal(False, vol(lambda_min-1.0))
        assert_equal(False, vol(lambda_max+1.0))
        # in state
        assert_equal(True, vol(lambda_min+1.0))
        assert_equal(True, vol(lambda_max-1.0))
        # border
        assert_equal(True, vol(lambda_min))
        assert_equal(True, vol(lambda_max))

    def test_thru_pbc_to_image(self):
        '''max in next periodic domain'''
        lambda_min = 70
        lambda_max = 210
        vol = volume.LambdaVolumePeriodic(op_id,
                                          lambda_min, lambda_max, -180,180)
        assert_equal(vol.lambda_max, -150)
        # assuming that's true, so is everything else

    @raises(Exception)
    def test_volume_bigger_than_bounds(self):
        '''max-min > pbc_range raises Exception'''
        vol = volume.LambdaVolumePeriodic(op_id, 90, 720, -180, 180)

    def test_volume_equals_bounds(self):
        '''max-min == pbc_range allows all points'''
        vol = volume.LambdaVolumePeriodic(op_id, 0, 360, -180, 180)
        assert_equal(vol.__str__(),
            "{x|(Id(x) - -180) % 360 + -180 in [-180, 180]}")
        assert_equal(True, vol(0))
        assert_equal(True, vol(360))
        assert_equal(True, vol(-180))
        assert_equal(True, vol(180))

