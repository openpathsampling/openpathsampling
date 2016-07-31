from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        assert_almost_equal, raises)
from nose.plugins.skip import Skip, SkipTest
from test_helpers import (
    true_func, assert_equal_array_array, make_1d_traj, data_filename
)

import openpathsampling as paths
from openpathsampling.high_level.interface_set import GenericVolumeInterfaceSet

import logging
logging.getLogger('openpathsampling.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.ensemble').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.netcdfplus').setLevel(logging.CRITICAL)

class testInterfaceSet(object):
    def setup(self):
        self.cv = paths.CV_Function(name="x", f=lambda s: s.xyz[0][0])
        self.lambdas = [0.0, 0.1, 0.2, 0.3]
        self.volumes = paths.VolumeFactory.CVRangeVolumeSet(self.cv, 
                                                            float("-inf"),
                                                            self.lambdas)
        self.interface_set = paths.InterfaceSet(self.volumes, self.cv,
                                                self.lambdas)
        self.no_lambda_set = paths.InterfaceSet(self.volumes, self.cv)

    def test_get_lambda(self):
        for (v, l) in zip(self.volumes, self.lambdas):
            assert_equal(self.interface_set.get_lambda(v), l)
            assert_equal(self.no_lambda_set.get_lambda(v), None)

    def test_list_behavior(self):
        # len
        assert_equal(len(self.interface_set), 4)
        assert_equal(len(self.no_lambda_set), 4)
        # getitem, contains
        for i in range(4):
            assert_equal(self.volumes[i], self.interface_set[i])
            assert_equal(self.volumes[i] in self.interface_set, True)
        # special case of -1 needs to work (used frequently!)
        assert_equal(self.volumes[-1], self.interface_set[-1])
        # iter
        for vol in self.interface_set:
            assert_equal(vol in self.volumes, True)
        # reversed
        i = 0 
        for vol in reversed(self.interface_set):
            assert_equal(vol, self.volumes[3-i])
            i += 1


class testGenericVolumeInterfaceSet(object):
    def test_determine_direction(self):
        # this is just to make the rest a little more readable
        determine_direction = GenericVolumeInterfaceSet._determine_direction
        assert_equal(1, determine_direction(float("-inf"), [0.0, 0.1, 0.2]))
        assert_equal(-1, determine_direction([0.2, 0.1, 0.0], float("inf")))
        assert_equal(0, determine_direction([-0.1, -0.2], [0.1, 0.2]))
        assert_equal(1, determine_direction([0.0, 0.0], [0.1, 0.2]))
        assert_equal(-1, determine_direction([-0.1, -0.2], [0.0, 0.0]))
        # and the idiot case:
        assert_equal(0, determine_direction([-0.1, -0.1], [0.1, 0.1]))

    @raises(RuntimeError)
    def test_bad_determine_direction(self):
        GenericVolumeInterfaceSet._determine_direction([0.0, -0.1], 
                                                       [0.1, 0.2, 0.3])
        

    def test_prep_minvals_maxvals(self):
        raise SkipTest


class testVolumeInterfaceSet(object):
    def setup(self):
        pass

    def test_get_lambda(self):
        raise SkipTest

    def test_list_behavior(self):
        raise SkipTest

    def test_new_interface(self):
        raise SkipTest
    

class testPeriodicVolumeInterfaceSet(object):
    def setup(self):
        pass

    def test_get_lambda(self):
        raise SkipTest

    def test_list_behavior(self):
        raise SkipTest

    def test_new_interface(self):
        raise SkipTest
    
