from nose.tools import assert_equal, assert_not_equal, raises
from nose.plugins.skip import SkipTest
from test_helpers import CallIdentity, prepend_exception_message, make_1d_traj


import openpathsampling as paths
from openpathsampling.ensemble import *

import logging
logging.getLogger('opentis.analysis.tis_analysis').setLevel(logging.DEBUG)
logging.getLogger('opentis.initialization').setLevel(logging.CRITICAL)

class testTISTransition(object):
    def setup(self):
        pass

    def test_initialization(self):
        pass

    def test_ensemble_statistics(self):
        pass

class testSummarizeTrajectoryVolumes(object):
    def setup(self):
        op = paths.CV_Function("Id", lambda snap : snap.coordinates[0][0])
        vol1 = paths.LambdaVolume(op, 0.1, 0.5)
        vol2 = paths.LambdaVolume(op, -0.1, 0.7)
        vol3 = paths.LambdaVolume(op, 2.0, 2.5)

        self.stateA = vol1
        self.interstitial = vol2 & ~vol1
        self.outInterface = ~vol2 & ~vol3
        self.stateB = vol3


    def _make_traj(self, string):
        pretraj = []
        for l in string:
            pretraj.append({"a" : 0.3, "b" : 2.2, "i" : 0.6, "x" : 1.0}[l])
        return make_1d_traj(coordinates=pretraj, velocities=[1.0]*len(pretraj))

    def test_summary_trajectory_volumes(self):
        pass

