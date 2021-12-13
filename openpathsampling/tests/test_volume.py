"""
@author David W.H. Swenson
"""
from __future__ import absolute_import

from builtins import object
from .test_helpers import (CallIdentity, raises_with_message_like,
                           make_1d_traj)

import pytest
import numpy as np

from openpathsampling.integration_tools import unit, HAS_SIMTK_UNIT

import openpathsampling.volume as volume

import openpathsampling as paths


class Identity2(CallIdentity):
    def __str__(self):
        return "Id2"


def setup_module():
    global op_id, volA, volB, volC, volD, volA2
    op_id = CallIdentity()
    volA = volume.CVDefinedVolume(op_id, -0.5, 0.5)
    volB = volume.CVDefinedVolume(op_id, 0.25, 0.75)
    volC = volume.CVDefinedVolume(op_id, -0.75, -0.25)
    volD = volume.CVDefinedVolume(op_id, -0.75, 0.75)
    volA2 = volume.CVDefinedVolume(Identity2(), -0.5, 0.5)


class TestEmptyVolume(object):
    def test_empty_volume(self):
        """Empty volume is well-behaved"""
        empty = volume.EmptyVolume()
        test = 0.1
        assert empty(test) is False
        assert (empty & volA)(test) is False
        assert (volA & empty)(test) is False
        assert (empty | volA)(test) is True
        assert (volA | empty)(test) is True
        assert (empty & volA).__str__() == "empty"
        # assert_is: logical combos with empty should return same obj
        assert (empty - volA) is empty
        assert (volA - empty) is volA
        assert (volA | empty) is volA
        assert (empty | volA) is volA
        assert (volA & empty) is empty
        assert (empty & volA) is empty
        assert (empty ^ volA) is volA
        assert (volA ^ empty) is volA
        assert (~ empty).__str__() == "all"

    def test_empty_volume_equality(self):
        empty1 = volume.EmptyVolume()
        empty2 = volume.EmptyVolume()
        assert empty1 == empty2
        assert not (empty1 != empty2)


class TestFullVolume(object):
    def test_full_volume(self):
        """Full volume is well-behaved"""
        full = volume.FullVolume()
        test = 0.1
        assert full(test) is True
        assert (full & volA)(test) is True
        assert (volA & full)(test) is True
        assert (full | volA)(test) is True
        assert (volA | full)(test) is True
        # assert_is: logical combos with full should return same obj
        assert (full & volA) is volA
        assert (volA & full) is volA
        assert (full | volA) is full
        assert (volA | full) is full
        assert (volA - full) == volume.EmptyVolume()
        assert (full - volA) == ~ volA
        assert (full ^ volA) == ~ volA
        assert (volA ^ full) == ~ volA
        assert (volA | full).__str__() == "all"
        assert (~ full).__str__() == "empty"


class TestCVDefinedVolume(object):
    def test_upper_boundary(self):
        assert volA(0.49) is True
        assert volA(0.50) is False
        assert volA(0.51) is False

    def test_inf_upper_boundary(self):
        volX = volume.CVDefinedVolume(op_id, 0, float('inf'))
        assert volA(float('inf')) is False
        assert volX(float('inf')) is True

    def test_lower_boundary(self):
        assert volA(-0.49) is True
        assert volA(-0.50) is True
        assert volA(-0.51) is False

    def test_inf_lower_boundary(self):
        volX = volume.CVDefinedVolume(op_id, float('-inf'), 0)
        assert volA(float('-inf')) is False
        assert volX(float('-inf')) is True

    def test_negation(self):
        assert (~volA)(0.25) is False
        assert (~volA)(0.75) is True
        assert (~volA)(0.5) is True
        assert (~volA)(-0.5) is False

    def test_autocombinations(self):
        # volA tests this in the CVRangeVolumes
        assert volA | volA is volA
        assert volA & volA is volA
        assert volA ^ volA == volume.EmptyVolume()
        assert volA - volA == volume.EmptyVolume()
        # combo tests this with VolumeCombination of CVRangeVolumes
        combo = (volD - volA)
        assert combo | combo is combo
        assert combo & combo is combo
        assert combo ^ combo == volume.EmptyVolume()
        assert combo - combo == volume.EmptyVolume()

    def test_and_combinations(self):
        assert (volA & volB) == volume.CVDefinedVolume(op_id, 0.25, 0.5)
        assert (volA & volB)(0.45) is True
        assert (volA & volB)(0.55) is False
        assert (volB & volC) == volume.EmptyVolume()
        # go to VolumeCombination if order parameter isn't the same
        assert (volA & volA2) == volume.VolumeCombination(volA, volA2,
                                                          lambda a, b: a and b,
                                                          '{0} and {1}')

    def test_or_combinations(self):
        assert (volA | volB) == volume.CVDefinedVolume(op_id, -0.5, 0.75)
        assert (volB | volC) == volume.UnionVolume(volB, volC)
        assert (volB | volC)(0.0) is False
        assert (volB | volC)(0.5) is True
        assert (volB | volC)(-0.5) is True

        # go to VolumeCombination if order parameters isn't the same
        assert (volA2 | volB) == volume.UnionVolume(volA2, volB)

    def test_xor_combinations(self):
        assert ((volA ^ volB) ==
                volume.UnionVolume(volume.CVDefinedVolume(op_id, -0.5, 0.25),
                                   volume.CVDefinedVolume(op_id, 0.5, 0.75)))
        assert (volA ^ volA2) == volume.SymmetricDifferenceVolume(volA, volA2)

    def test_sub_combinations(self):
        assert (volA - volB) == volume.CVDefinedVolume(op_id, -0.5, 0.25)
        assert (volB - volC) == volB
        assert (volA - volD) == volume.EmptyVolume()
        assert (volB - volA) == volume.CVDefinedVolume(op_id, 0.5, 0.75)
        assert ((volD - volA) ==
                volume.UnionVolume(volume.CVDefinedVolume(op_id, -0.75, -0.5),
                                   volume.CVDefinedVolume(op_id, 0.5, 0.75)))
        assert (volA2 - volA) == volume.RelativeComplementVolume(volA2, volA)

    def test_str(self):
        assert volA.__str__() == "{x|Id(x) in [-0.5, 0.5]}"
        assert (~volA).__str__() == "(not {x|Id(x) in [-0.5, 0.5]})"

    @pytest.mark.skipif(not paths.integration_tools.HAS_SIMTK_UNIT,
                        reason="Needs openmm.unit")
    def test_unit_support(self):
        u = unit

        vol = volume.CVDefinedVolume(
            op_id, -0.5 * u.nanometers, 0.25 * u.nanometers)

        assert vol(-0.25 * u.nanometers)
        assert not vol(-0.75 * u.nanometers)

        vol = volume.PeriodicCVDefinedVolume(
            op_id,
            -30 * u.nanometers, 90 * u.nanometers,
            -180 * u.nanometers, 180 * u.nanometers)

        assert vol(50 * u.nanometers)
        assert not vol(-70 * u.nanometers)

    @staticmethod
    def _vol_for_cv_type(inp):
        if not HAS_SIMTK_UNIT and inp == 'simtk':
            pytest.skip()

        func = {
            'float': lambda s: 1.0,
            'array': lambda s: np.array([1.0, 2.0]),
            'array1': lambda s: np.array([1.0]),
            'simtk': None
        }[inp]
        if func is None:  # only if inp is 'simtk'
            func = lambda s: 1.0 * unit.nanometers

        cv = paths.FunctionCV('cv', func)
        volume = paths.CVDefinedVolume(cv, 0.0, 1.5)
        return volume

    @pytest.mark.parametrize('inp', ['float', 'array', 'array1', 'simtk'])
    def test_is_iterable(self, inp):
        snap = make_1d_traj([0.0])[0]
        volume = self._vol_for_cv_type(inp)
        val = volume.collectivevariable(snap)
        expected = inp in ['array', 'array1']
        if expected:
            with pytest.warns(UserWarning, match="returns an iterable"):
                result = volume._is_iterable(val)
        else:
            result = volume._is_iterable(val)

        assert result is expected

    @pytest.mark.parametrize('inp', ['float', 'array1', 'simtk'])
    @pytest.mark.filterwarnings("ignore:The CV 'cv' returns an iterable")
    def test_get_cv_float(self, inp):
        snap = make_1d_traj([0.0])[0]
        volume = self._vol_for_cv_type(inp)
        val = volume._get_cv_float(snap)
        expected = inp in ['float', 'array1']
        assert isinstance(val, float) is expected


class TestCVRangeVolumePeriodic(object):
    def setup(self):
        self.pvolA = volume.PeriodicCVDefinedVolume(op_id, -100, 75)
        self.pvolA_ = volume.PeriodicCVDefinedVolume(op_id, 75, -100)
        self.pvolB = volume.PeriodicCVDefinedVolume(op_id, 50, 100)
        self.pvolC = volume.PeriodicCVDefinedVolume(op_id, -100, -50)
        self.pvolD = volume.PeriodicCVDefinedVolume(op_id, -100, 100)
        self.pvolE = volume.PeriodicCVDefinedVolume(op_id, -150, 150)

    def test_normal(self):
        """min<max and both within periodic domain"""
        lambda_min = -150
        lambda_max = 70
        vol = volume.PeriodicCVDefinedVolume(op_id, lambda_min, lambda_max,
                                             -180, 180)
        assert vol._period_len == 360
        assert vol._period_shift == -180
        assert (vol.__str__() ==
                "{x|(Id(x) - -180) % 360 + -180 in [-150, 70]}")
        for periodic_image in [-1, 0, 1]:
            # out of state
            assert vol(lambda_min-1.0 + 360*periodic_image) is False
            assert vol(lambda_max+1.0 + 360*periodic_image) is False
            # in state
            assert vol(lambda_min+1.0 + 360*periodic_image) is True
            assert vol(lambda_max-1.0 + 360*periodic_image) is True
            # border
            assert vol(lambda_min + 360*periodic_image) is True
            assert vol(lambda_max + 360*periodic_image) is False

    def test_normal_implicit(self):
        """min<max, no periodic domain defined"""
        lambda_min = -150
        lambda_max = 70
        vol = volume.PeriodicCVDefinedVolume(op_id,
                                             lambda_min, lambda_max)
        assert vol.__str__() == "{x|Id(x) [periodic] in [-150, 70]}"
        # out of state
        assert vol(lambda_min-1.0) is False
        assert vol(lambda_max+1.0) is False
        # in state
        assert vol(lambda_min+1.0) is True
        assert vol(lambda_max-1.0) is True
        # border
        assert vol(lambda_min) is True
        assert vol(lambda_max) is False

    def test_wrapped_volume(self):
        """max<min and both within periodic domain"""
        lambda_min = 70
        lambda_max = -150
        vol = volume.PeriodicCVDefinedVolume(op_id,
                                             lambda_min, lambda_max, -180, 180)
        assert vol.__str__() == ("{x|(Id(x) - -180) % 360 + -180 "
                                 "in [-180, -150] union [70, 180]}")
        for periodic_image in [-1, 0, 1]:
            # out of state
            assert vol(lambda_min-1.0 + 360*periodic_image) is False
            assert vol(lambda_max+1.0 + 360*periodic_image) is False
            # in state
            assert vol(lambda_min+1.0 + 360*periodic_image) is True
            assert vol(lambda_max-1.0 + 360*periodic_image) is True
            # border
            assert vol(lambda_min + 360*periodic_image) is True
            assert vol(lambda_max + 360*periodic_image) is False

    def test_wrapped_volume_implicit(self):
        """max<min, no periodic domain defined"""
        lambda_min = 70
        lambda_max = -150
        vol = volume.PeriodicCVDefinedVolume(op_id, lambda_min, lambda_max)
        # out of state
        assert vol(lambda_min-1.0) is False
        assert vol(lambda_max+1.0) is False
        # in state
        assert vol(lambda_min+1.0) is True
        assert vol(lambda_max-1.0) is True
        # border
        assert vol(lambda_min) is True
        assert vol(lambda_max) is False

    def test_thru_pbc_to_image(self):
        '''max in next periodic domain'''
        lambda_min = 70
        lambda_max = 210
        vol = volume.PeriodicCVDefinedVolume(op_id, lambda_min, lambda_max,
                                             -180, 180)
        assert vol.lambda_max == -150
        # assuming that's true, so is everything else

    def test_volume_bigger_than_bounds(self):
        '''max-min > pbc_range raises Exception'''
        with pytest.raises(Exception, match='Range of volume'):
            volume.PeriodicCVDefinedVolume(op_id, 90, 720, -180, 180)

    def test_volume_equals_bounds(self):
        '''max-min == pbc_range allows all points'''
        vol = volume.PeriodicCVDefinedVolume(op_id, 0, 360, -180, 180)
        assert vol.__str__() == ("{x|(Id(x) - -180) % 360 + -180 "
                                 "in [-180, 180]}")
        assert vol(0) is True
        assert vol(360) is True
        assert vol(-180) is True
        assert vol(180) is True

    def test_periodic_and_combos(self):
        assert ((self.pvolA & self.pvolB) ==
                volume.PeriodicCVDefinedVolume(op_id, 50, 75))
        assert (self.pvolA & self.pvolB)(60) is True
        assert (self.pvolA & self.pvolB)(80) is False
        assert (self.pvolB & self.pvolC) == volume.EmptyVolume()
        assert (self.pvolC & self.pvolB) == volume.EmptyVolume()
        assert (self.pvolA & self.pvolA) is self.pvolA
        assert (self.pvolA & self.pvolA_) == volume.EmptyVolume()
        assert (self.pvolE & self.pvolD) == self.pvolD
        # go to special case for cyclic permutation
        assert (self.pvolB & self.pvolD) == self.pvolB
        # go to special case
        assert ((self.pvolE & self.pvolA_) ==
                volume.UnionVolume(
                    volume.PeriodicCVDefinedVolume(op_id, -150, -100),
                    volume.PeriodicCVDefinedVolume(op_id, 75, 150)
                    )
                )
        # go to super if needed
        assert type(self.pvolA & volA) is volume.IntersectionVolume

    def test_periodic_or_combos(self):
        assert (self.pvolA | self.pvolB) == self.pvolD
        assert (self.pvolA | self.pvolB)(60) is True
        assert (self.pvolA | self.pvolB)(80) is True
        assert (self.pvolA | self.pvolB)(125) is False
        assert ((self.pvolB | self.pvolC) ==
                volume.UnionVolume(self.pvolB, self.pvolC))
        assert ((self.pvolC | self.pvolB) ==
                volume.UnionVolume(self.pvolC, self.pvolB))
        assert (self.pvolA | self.pvolA) is self.pvolA
        assert (self.pvolA | self.pvolA_) == volume.FullVolume()
        assert (self.pvolE | self.pvolD) == self.pvolE
        assert (self.pvolB | self.pvolD) == self.pvolD
        assert (self.pvolE | self.pvolA_) == volume.FullVolume()

    def test_periodic_xor_combos(self):
        assert self.pvolA ^ self.pvolA_ == volume.FullVolume()
        assert self.pvolA ^ self.pvolA == volume.EmptyVolume()
        assert (self.pvolE ^ self.pvolD ==
                volume.UnionVolume(
                    volume.PeriodicCVDefinedVolume(op_id, -150, -100),
                    volume.PeriodicCVDefinedVolume(op_id, 100, 150)))
        assert self.pvolB ^ self.pvolC == self.pvolB | self.pvolC
        assert (self.pvolB ^ self.pvolD ==
                volume.PeriodicCVDefinedVolume(op_id, -100, 50))

    def test_periodic_not_combos(self):
        assert ~self.pvolA == self.pvolA_
        assert self.pvolA == ~self.pvolA_
        assert (~(self.pvolB | self.pvolC) ==
                volume.NegatedVolume(self.pvolB | self.pvolC))
        assert (~(self.pvolB | self.pvolC))(25) is True
        assert (~(self.pvolB | self.pvolC))(75) is False
        assert (~(self.pvolB | self.pvolC))(-75) is False

    def test_periodic_sub_combos(self):
        assert self.pvolA - self.pvolA_ is self.pvolA
        assert self.pvolA_ - self.pvolA is self.pvolA_
        assert (self.pvolD - self.pvolA ==
                volume.PeriodicCVDefinedVolume(op_id, 75, 100))
        assert (self.pvolD - self.pvolA)(80) is True
        assert (self.pvolD - self.pvolA)(50) is False
        assert self.pvolB - self.pvolC is self.pvolB
        assert self.pvolA - self.pvolA == volume.EmptyVolume()
        assert (self.pvolE - self.pvolD ==
                volume.UnionVolume(
                    volume.PeriodicCVDefinedVolume(op_id, -150, -100),
                    volume.PeriodicCVDefinedVolume(op_id, 100, 150)))
        assert (self.pvolE - self.pvolA_ ==
                volume.PeriodicCVDefinedVolume(op_id, -100, 75))


class TestAbstract(object):
    @raises_with_message_like(TypeError, "Can't instantiate abstract class")
    def test_abstract_volume(self):
        volume.Volume()
