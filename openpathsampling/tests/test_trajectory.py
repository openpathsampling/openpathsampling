from __future__ import absolute_import
from builtins import object
import logging

from nose.tools import (
    assert_equal, assert_not_equal, raises
)
from nose.plugins.skip import SkipTest
from .test_helpers import (CallIdentity, prepend_exception_message,
                           make_1d_traj, assert_items_equal)


import openpathsampling as paths
from .test_helpers import make_1d_traj

logging.getLogger('opentis.trajectory').setLevel(logging.DEBUG)
logging.getLogger('opentis.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)

class TestSummarizeTrajectoryVolumes(object):
    def setup(self):
        op = paths.FunctionCV("Id", lambda snap : snap.coordinates[0][0])
        vol1 = paths.CVDefinedVolume(op, 0.1, 0.5)
        vol2 = paths.CVDefinedVolume(op, -0.1, 0.7)
        vol3 = paths.CVDefinedVolume(op, 2.0, 2.5)

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

    def test_summarize_trajectory_volumes_str(self):
        voldict = {"A" : self.stateA, "B" : self.stateB,
                   "I" : self.interstitial, "X" : self.outInterface}
        assert_items_equal(
            self._make_traj("abix").summarize_by_volumes_str(voldict),
            "A-B-I-X"
        )
        assert_items_equal(
            self._make_traj("abix").summarize_by_volumes_str(voldict, " blah "), 
            "A blah B blah I blah X"
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

    def test_summarize_trajectory_volumes_with_nonevol(self):
        voldict = {"A" : self.stateA, "B" : self.stateB,
                   "I" : self.interstitial}
        assert_items_equal(
            self._make_traj("abix").summarize_by_volumes_str(voldict),
            "A-B-I-None"
        )
        assert_items_equal(
            self._make_traj("aiiibbxxbx").summarize_by_volumes_str(voldict),
            "A-I-B-None-B-None"
        )

class TestSubtrajectoryIndices(object):
    def setup(self):
        op = paths.FunctionCV("Id", lambda snap : snap.coordinates[0][0])
        vol1 = paths.CVDefinedVolume(op, 0.1, 0.5)
        vol2 = paths.CVDefinedVolume(op, -0.1, 0.7)
        vol3 = paths.CVDefinedVolume(op, 2.0, 2.5)

        self.stateA = vol1
        self.interstitial = vol2 & ~vol1
        self.outInterface = ~vol2 & ~vol3
        self.stateB = vol3

    def test_subtrajectory_indices(self):
        # simplify more complicated expressions
        stateA = self.stateA
        stateB = self.stateB
        pretraj = [0.20, 0.30, 0.60, 0.40, 0.65, 2.10, 2.20, 2.60, 2.10,
                   0.80, 0.55, 0.40, 0.20]
        # 00, 01, 02, 03, 04, 05, 06, 07, 08, 09, 10, 11, 12
        #  A,  A,  I,  A,  I,  B,  B,  X,  B,  X,  I,  A,  A
        trajectory = make_1d_traj(
            coordinates=pretraj,
            velocities=[1.0]*len(pretraj)
        )
        ensemble_A = paths.AllInXEnsemble(stateA)
        ensemble_B = paths.AllInXEnsemble(stateB)
        ensemble_ABA = paths.SequentialEnsemble([
            paths.AllInXEnsemble(stateA) & paths.LengthEnsemble(1),
            paths.PartInXEnsemble(stateB) & paths.AllOutXEnsemble(stateA),
            paths.AllInXEnsemble(stateA) & paths.LengthEnsemble(1)
        ])
        subtrajectoriesA = ensemble_A.split(trajectory, overlap=0)
        subtrajectoriesB = ensemble_B.split(trajectory, overlap=0)
        subtrajectoriesABA = ensemble_ABA.split(trajectory)

        # make sure we have the trajectories we expect
        assert_equal(len(subtrajectoriesA), 3)
        assert_equal(len(subtrajectoriesB), 2)
        assert_equal(len(subtrajectoriesABA), 1)
        # the following assertions check that the subtrajectories are the
        # ones that we expect; the numbers here are linked to the indices
        # we'll test next
        assert_equal(subtrajectoriesA[0], trajectory[0:2])
        assert_equal(subtrajectoriesA[1], trajectory[3:4])
        assert_equal(subtrajectoriesA[2], trajectory[11:13])
        assert_equal(subtrajectoriesB[0], trajectory[5:7])
        assert_equal(subtrajectoriesB[1], trajectory[8:9])
        assert_equal(subtrajectoriesABA[0], trajectory[3:12])
        # now we run the subtrajectory_indices function and test it
        indicesA = trajectory.subtrajectory_indices(subtrajectoriesA)
        indicesB = trajectory.subtrajectory_indices(subtrajectoriesB)
        indicesABA = trajectory.subtrajectory_indices(subtrajectoriesABA)
        assert_equal(indicesA, [[0, 1], [3], [11, 12]])
        assert_equal(indicesB, [[5, 6], [8]])
        assert_equal(indicesABA, [[3, 4, 5, 6, 7, 8, 9, 10, 11]])
