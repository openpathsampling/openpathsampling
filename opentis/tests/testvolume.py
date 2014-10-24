"""
@author David W.H. Swenson
"""

import os
import sys
from nose.tools import assert_equal, assert_not_equal, assert_is, raises
from nose.plugins.skip import Skip, SkipTest
from test_helpers import CallIdentity

sys.path.append(os.path.abspath('../'))
import volume

class Identity2(CallIdentity):
    def __str__(self):
        return "Id2"

def setUp():
    global op_id, volA, volB, volC, volD, volA2
    op_id = CallIdentity()
    volA = volume.LambdaVolume(op_id, -0.5, 0.5)
    volB = volume.LambdaVolume(op_id, 0.25, 0.75)
    volC = volume.LambdaVolume(op_id, -0.75, -0.25)
    volD = volume.LambdaVolume(op_id, -0.75, 0.75)
    volA2 = volume.LambdaVolume(Identity2(), -0.5, 0.5)


class testEmptyVolume(object):
    def test_empty_volume(self):
        """Empty volume is well-behaved"""
        empty = volume.EmptyVolume()
        test = 0.1
        assert_equal(empty(test), False)
        assert_equal((empty & volA)(test), False)
        assert_equal((volA & empty)(test), False)
        assert_equal((empty | volA)(test), True)
        assert_equal((volA | empty)(test), True)
        assert_equal((empty & volA).__str__(), "empty")
        # assert_is: logical combos with empty should return same obj
        assert_is((empty - volA), empty)
        assert_is((volA - empty), volA)
        assert_is((volA | empty), volA)
        assert_is((empty | volA), volA)
        assert_is((volA & empty), empty)
        assert_is((empty & volA), empty)
        assert_is((empty ^ volA), volA)
        assert_is((volA ^ empty), volA)
        assert_equal((~ empty).__str__(), "all")

class testFullVolume(object):
    def test_full_volume(self):
        """Full volume is well-behaved"""
        full = volume.FullVolume()
        test = 0.1
        assert_equal(full(test), True)
        assert_equal((full & volA)(test), True)
        assert_equal((volA & full)(test), True)
        assert_equal((full | volA)(test), True)
        assert_equal((volA | full)(test), True)
        # assert_is: logical combos with full should return same obj
        assert_is((full & volA), volA)
        assert_is((volA & full), volA)
        assert_is((full | volA), full)
        assert_is((volA | full), full)
        assert_equal((volA - full), volume.EmptyVolume())
        assert_equal((full - volA), ~ volA)
        assert_equal((full ^ volA), ~ volA)
        assert_equal((volA ^ full), ~ volA)
        assert_equal((volA | full).__str__(), "all")
        assert_equal((~ full).__str__(), "empty")

class testLambdaVolume(object):
    def setUp(self):
        pass

    def teardown(self):
        pass

    def test_lower_boundary(self):
        pass

    def test_upper_boundary(self):
        pass

    def test_negation(self):
        pass

    def test_autocombinations(self):
        assert_equal(volA | volA, volA)
        assert_equal(volA & volA, volA)
        assert_equal(volA ^ volA, volume.EmptyVolume())
        assert_equal(volA - volA, volume.EmptyVolume())

    def test_and_combinations(self):
        assert_equal((volA & volB), volume.LambdaVolume(op_id, 0.25, 0.5))
        assert_equal((volA & volB)(0.45), True)
        assert_equal((volA & volB)(0.55), False)
        assert_equal((volB & volC), volume.EmptyVolume())
        # go to VolumeCombination if order parameter isn't the same
        assert_equal((volA & volA2),
                     volume.VolumeCombination(volA, volA2,
                                              lambda a, b: a and b,
                                              '{0} and {1}')
                    )

    def test_or_combinations(self):
        or_fnc = lambda a, b : a or b
        or_str = '{0} or {1}'
        assert_equal((volA | volB), volume.LambdaVolume(op_id, -0.5, 0.75))
        assert_equal((volB | volC), 
                     volume.VolumeCombination(volB, volC, or_fnc, or_str))
        assert_equal((volB | volC)(0.0), False)
        assert_equal((volB | volC)(0.5), True)
        assert_equal((volB | volC)(-0.5), True)
        
        # go to VolumeCombination if order parameters isn't the same
        assert_equal((volA2 | volB),
                     volume.VolumeCombination(volA2, volB, or_fnc, or_str))

    def test_xor_combinations(self):
        or_fcn = lambda a, b: a | b
        or_str = '{0} or {1}'
        xor_fcn = lambda a, b: a ^ b
        xor_str = '{0} xor {1}'

        assert_equal((volA ^ volB),
                     volume.VolumeCombination(
                         volume.LambdaVolume(op_id, -0.5, 0.25),
                         volume.LambdaVolume(op_id, 0.5, 0.75),
                         or_fcn, or_str
                     ))
        assert_equal((volA ^ volA2),
                     volume.VolumeCombination(volA, volA2, xor_fcn, xor_str))

        pass

    def test_sub_combinations(self):
        sub_fcn = lambda a, b: a and not b
        sub_str = '{0} and not {1}'
        assert_equal((volA - volB), volume.LambdaVolume(op_id, -0.5, 0.25))
        assert_equal((volB - volC), volB)
        assert_equal((volA - volD), volume.EmptyVolume())
        assert_equal((volB - volA), volume.LambdaVolume(op_id, 0.5, 0.75))
        assert_equal((volD - volA),
                     volume.VolumeCombination(
                         volume.LambdaVolume(op_id, -0.75, -0.5),
                         volume.LambdaVolume(op_id, 0.5, 0.75),
                         lambda a, b: a or b, '{0} or {1}')
                    )
        assert_equal((volA2 - volA),
                     volume.VolumeCombination(volA2, volA, sub_fcn, sub_str))
        pass

    def test_str(self):
        # string and inverted string
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

