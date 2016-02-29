from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        raises, assert_almost_equal)
from nose.plugins.skip import SkipTest
from test_helpers import make_1d_traj

import openpathsampling as paths

import random

import logging
logging.getLogger('openpathsampling.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.netcdfplus').setLevel(logging.CRITICAL)

class testSingleTrajectoryAnalysis(object):
    def setup(self):
        op = paths.CV_Function("Id", lambda snap : snap.coordinates[0][0])
        vol1 = paths.CVRangeVolume(op, 0.1, 0.5)
        vol3 = paths.CVRangeVolume(op, 2.0, 2.5)
        self.stateA = vol1
        self.stateB = vol3

        self.stateX = ~vol1 & ~vol3

        transition = paths.TPSTransition(self.stateA, self.stateB)
        self.analyzer = paths.SingleTrajectoryAnalysis(transition, dt=0.1)
        self.traj_str = "aaaxaxxbxaxababaxbbbbbxxxxxxa"
        # frame numbers "0    5    0    5    0    5  8"
        self.trajectory = self._make_traj(self.traj_str)

    def _make_traj(self, traj_str):
        sequence = []
        char_to_parameters = {
            'a' : {'lo' : 0.1, 'hi' : 0.5},
            'b' : {'lo' : 2.0, 'hi' : 2.5},
            'x' : {'lo' : 0.5, 'hi' : 2.0}
        }
        delta = 0.05
        for char in traj_str:
            params = char_to_parameters[char]
            n_max = int((params['hi'] - params['lo'])/delta)
            sequence.append(params['lo'] + delta*random.randint(1, n_max-1))

        return make_1d_traj(coordinates=sequence,
                            velocities=[1.0]*len(sequence))

    def test_analyze_continuous_time(self):
        self.analyzer.analyze_continuous_time(self.trajectory, self.stateA)
        assert_equal(self.analyzer.continuous_frames[self.stateA].tolist(),
                     [3, 1, 1, 1, 1, 1, 1])
        self.analyzer.analyze_continuous_time(self.trajectory, self.stateB)
        assert_equal(self.analyzer.continuous_frames[self.stateB].tolist(),
                     [1, 1, 1, 5])
        assert_almost_equal(self.analyzer.continuous_times[self.stateA].mean(),
                            9.0/7.0*0.1)
        assert_almost_equal(self.analyzer.continuous_times[self.stateB].mean(),
                            8.0/4.0*0.1)

    @raises(KeyError)
    def test_analyze_continuous_time_bad_state(self):
        self.analyzer.analyze_continuous_time(self.trajectory, self.stateX)

    def test_analyze_lifetime(self):
        self.analyzer.analyze_lifetime(self.trajectory, self.stateA)
        self.analyzer.analyze_lifetime(self.trajectory, self.stateB)
        assert_equal(self.analyzer.lifetime_frames[self.stateA].tolist(),
                     [3, 1, 2])  # A->B
        assert_equal(self.analyzer.lifetime_frames[self.stateB].tolist(),
                     [2, 1, 1, 11])  # B->A
        assert_almost_equal(self.analyzer.lifetimes[self.stateA].mean(),
                            6.0/3.0*0.1)
        assert_almost_equal(self.analyzer.lifetimes[self.stateB].mean(),
                            15.0/4.0*0.1)

    def test_analyze_flux(self):
        # A: [{out: 1, in: 1}
        # B: insufficient 
        # NOTE: we may want a separate trajectory for this
        raise SkipTest

    def test_analyze(self):
        raise SkipTest

    def test_summary(self):
        raise SkipTest
