from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        assert_almost_equal, raises)

from nose.plugins.skip import Skip, SkipTest
from test_helpers import true_func, assert_equal_array_array, make_1d_traj

import openpathsampling as paths
from openpathsampling.bias_function import *

import logging
from openpathsampling import VolumeFactory as vf

logging.getLogger('openpathsampling.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.ensemble').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)

class testBiasEnsembleTable(object):
    def setup(self):
        # create the network
        xval = paths.CV_Function(name="xA", f=lambda s : s.xyz[0][0])
        self.stateA = paths.CVRangeVolume(xval, float("-inf"), -0.5)
        ifacesA = vf.CVRangeVolumeSet(xval, float("-inf"), [-0.5, -0.4, -0.3])
        self.network = paths.MISTISNetwork([
            (self.stateA, ifacesA, xval, self.stateA)
        ])
        # create the 
        pass

    def test_bias_ensemble_old_to_new(self):
        pass

    def test_bias_ensemble_new_to_old(self):
        pass


