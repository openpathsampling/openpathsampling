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

class testMSOuterTISInterface(object):
    def setup(self):
        self.cv = paths.CV_Function(name="x", f=lambda s: s.xyz[0][0])
        self.lambdas = [0.0, 0.1, 0.2, 0.3]
        self.interfaces = paths.VolumeInterfaceSet(cv=self.cv,
                                                   minvals=float("-inf"),
                                                   maxvals=self.lambdas)
        self.stateA = paths.CVRangeVolume(self.cv, float("-inf"), 0.0
                                         ).named("A")
        self.stateB = paths.CVRangeVolume(self.cv, 1.0, float("inf")
                                         ).named("B")
        self.network = paths.MISTISNetwork([
            (self.stateA, self.interfaces, self.stateB)
        ])
            

        # self.network
        # self.ms_outer
        pass

    def test_volume_for_interface_set(self):
        raise SkipTest

    def test_lambda_for_interface_set(self):
        raise SkipTest

    def test_make_ensemble(self):
        raise SkipTest
