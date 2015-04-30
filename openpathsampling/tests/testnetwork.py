import os
import numpy as np

from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        assert_almost_equal, raises)
from nose.plugins.skip import Skip, SkipTest
from test_helpers import true_func, assert_equal_array_array, make_1d_traj

import openpathsampling as paths
from openpathsampling.analysis.network import *
from openpathsampling import VolumeFactory as vf

class testMSTISNetwork(object):
    def setup(self):
        xval = paths.CV_Function(name="xA", fcn=lambda s : s.xyz[0][0])
        stateA = paths.LambdaVolume(xval, float("-inf"), -0.5)
        stateB = paths.LambdaVolume(xval, -0.1, 0.1)
        stateC = paths.LambdaVolume(xval, 0.5, float("inf"))

        ifacesA = vf.LambdaVolumeSet(xval, float("-inf"), [-0.5, -0.4, -0.3])
        ifacesB = vf.LambdaVolumeSet(xval, [-0.2, -0.15, -0.1], [0.2, 0.15, 0.1])
        ifacesC = vf.LambdaVolumeSet(xval, [0.5, 0.4, 0.3], float("inf"))

        self.mstis = MSTISNetwork([
            (stateA, ifacesA, "A", xval),
            (stateB, ifacesB, "B", xval),
            (stateC, ifacesC, "C", xval)
        ])


        self.traj = {}
        #self.traj['AA']
        #self.traj['AB']
        #self.traj['BA']
        #self.traj['BB']
        #self.traj['BC']
        #self.traj['CB']
        #self.traj['CC']
        # A->C magically jumps over B
        #self.traj['AC']
        #self.traj['CA']
        pass

    def test_trajectories(self):
        pass

    def test_ms_outers(self):
        pass
