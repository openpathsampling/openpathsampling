from nose.tools import assert_equal, assert_not_equal, assert_items_equal, raises
from nose.plugins.skip import SkipTest
from test_helpers import CallIdentity, prepend_exception_message, make_1d_traj


import openpathsampling as paths
from openpathsampling.trajectory import *

import logging
logging.getLogger('opentis.trajectory').setLevel(logging.DEBUG)
logging.getLogger('opentis.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)

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

    def test_summarize_trajectory_volumes(self):
        voldict = {"A" : self.stateA, "B" : self.stateB, 
                   "I" : self.interstitial, "X" : self.outInterface}
        assert_items_equal(
            self._make_traj("abix").summarize_by_volumes(voldict), 
            [("A", 1), ("B", 1), ("I", 1), ("X", 1)]
        )
        assert_items_equal(
            self._make_traj("aiiibbxxbx").summarize_by_volumes(voldict), 
            [("A", 1), ("I", 3), ("B", 2), ("X", 2), ("B", 1), ("X", 1)]
        )

    def test_summarize_trajectory_volumes_with_nonevol(self):
        voldict = {"A" : self.stateA, "B" : self.stateB, 
                   "I" : self.interstitial}
        assert_items_equal(
            self._make_traj("abix").summarize_by_volumes(voldict), 
            [("A", 1), ("B", 1), ("I", 1), (None, 1)]
        )
        assert_items_equal(
            self._make_traj("aiiibbxxbx").summarize_by_volumes(voldict), 
            [("A", 1), ("I", 3), ("B", 2), (None, 2), ("B", 1), (None, 1)]
        )


