from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        assert_almost_equal, raises)
from nose.plugins.skip import Skip, SkipTest
from test_helpers import true_func, assert_equal_array_array, make_1d_traj


import openpathsampling as paths
from openpathsampling.analysis.move_strategy import *
from openpathsampling import VolumeFactory as vf

import logging
logging.getLogger('openpathsampling.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.ensemble').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)

class testMoveStrategy(object):
    def setup(self):
        cvA = paths.CV_Function(name="xA", fcn=lambda s : s.xyz[0][0])
        cvB = paths.CV_Function(name="xA", fcn=lambda s : -s.xyz[0][0])
        stateA = paths.LambdaVolume(cvA, float("-inf"), -0.5)
        stateB = paths.LambdaVolume(cvB, float("-inf"), -0.5)
        interfacesA = vf.LambdaVolumeSet(cvA, float("-inf"), [-0.5, -0.3])
        interfacesB = vf.LambdaVolumeSet(cvB, float("-inf"), [-0.5, -0.3])
        self.network = paths.MSTISNetwork([
            (stateA, interfacesA, "A", cvA),
            (stateB, interfacesB, "B", cvB)
        ])

    def test_make_chooser(self):
        raise SkipTest

    def test_get_scheme(self):
        raise SkipTest

    def test_get_ensembles(self):
        raise SkipTest
