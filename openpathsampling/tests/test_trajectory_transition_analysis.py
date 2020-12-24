from __future__ import division
from __future__ import absolute_import
from builtins import zip
from past.utils import old_div
from builtins import object
from nose.tools import (assert_equal, assert_not_equal, raises,
                        assert_almost_equal)
from nose.plugins.skip import SkipTest
from .test_helpers import make_1d_traj

import openpathsampling as paths

import random
import numpy as np

import logging
logging.getLogger('openpathsampling.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.ensemble').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.netcdfplus').setLevel(logging.CRITICAL)

class TestTrajectorySegmentContainer(object):
    def setup(self):
        op = paths.FunctionCV("Id", lambda snap : snap.coordinates[0][0])
        self.vol1 = paths.CVDefinedVolume(op, 0.1, 0.5)
        self.vol3 = paths.CVDefinedVolume(op, 2.0, 2.5)

        self.trajectory = make_1d_traj(coordinates=[0.2, 0.3, 0.6, 2.1, 2.2,
                                                    0.7, 0.4, 0.35, 2.4,
                                                    0.33, 0.32, 0.31],
                                       velocities=[0.0]*12)

        all_in_1 = paths.AllInXEnsemble(self.vol1)
        self.segments = all_in_1.split(self.trajectory)
        self.container = paths.TrajectorySegmentContainer(self.segments,
                                                          dt=0.5)

    def test_list_behavior(self):
        assert_equal(len(self.container), 3)
        assert_equal(self.container[0], self.segments[0])
        for (truth, beauty) in zip(self.container, self.segments):
            assert_equal(truth, beauty)
        assert_equal(self.segments[0] in self.container, True)

    @raises(TypeError)
    def test_segments_setitem_fails(self):
        self.container[0] = self.trajectory

    def test_segments(self):
        assert_equal(self.container[0], self.trajectory[0:2])
        assert_equal(self.container[1], self.trajectory[6:8])
        assert_equal(self.container[2], self.trajectory[9:12])

    def test_from_trajectory_and_indices(self):
        container = \
            paths.TrajectorySegmentContainer.from_trajectory_and_indices(
                trajectory=self.trajectory,
                indices=[(0, 2), (6, 8), (9, 12)],
                dt=0.5
            )
        assert_equal(container[0], self.trajectory[0:2])
        assert_equal(container[1], self.trajectory[6:8])
        assert_equal(container[2], self.trajectory[9:12])

    def test_n_frames(self):
        assert_equal(self.container.n_frames.tolist(), [2, 2, 3])

    def test_times(self):
        assert_equal(self.container.times.tolist(), [1.0, 1.0, 1.5])

    @raises(RuntimeError)
    def test_times_without_dt(self):
        bad_container = paths.TrajectorySegmentContainer(self.segments)
        bad_container.times

    def test_add(self):
        ens_B = paths.AllInXEnsemble(self.vol3)
        segs_B = ens_B.split(self.trajectory)
        container_B = paths.TrajectorySegmentContainer(segs_B, dt=0.5)
        assert_equal(len(container_B), 2)
        test_container = self.container + container_B
        assert_equal(len(test_container), 5)
        for seg in self.container:
            assert_equal(seg in test_container, True)
        for seg in container_B:
            assert_equal(seg in test_container, True)

    @raises(RuntimeError)
    def test_add_different_dt(self):
        ens_B = paths.AllInXEnsemble(self.vol3)
        segs_B = ens_B.split(self.trajectory)
        container_B = paths.TrajectorySegmentContainer(segs_B)
        test_container = self.container + container_B

    def test_iadd(self):
        ens_B = paths.AllInXEnsemble(self.vol3)
        segs_B = ens_B.split(self.trajectory)
        container_B = paths.TrajectorySegmentContainer(segs_B)
        container_B_id = id(container_B)
        container_B += self.container
        assert_equal(len(container_B), 5)
        assert_equal(container_B_id, id(container_B))


class TestTrajectoryTransitionAnalysis(object):
    def setup(self):
        op = paths.FunctionCV("Id", lambda snap : snap.coordinates[0][0])
        vol1 = paths.CVDefinedVolume(op, 0.1, 0.5)
        vol2 = paths.CVDefinedVolume(op, -0.1, 0.7)
        vol3 = paths.CVDefinedVolume(op, 2.0, 2.5)
        self.stateA = vol1
        self.stateB = vol3
        self.interfaceA0 = vol2

        self.stateX = ~vol1 & ~vol3

        self.transition = paths.TPSTransition(self.stateA, self.stateB)
        self.analyzer = paths.TrajectoryTransitionAnalysis(self.transition,
                                                           dt=0.1)
        self.traj_str = "aaaxaxxbxaxababaxbbbbbxxxxxxa"
        # frame numbers "0    5    0    5    0    5  8"
        self.trajectory = self._make_traj(self.traj_str)

    def _make_traj(self, traj_str):
        sequence = []
        char_to_parameters = {
            'a' : {'lo' : 0.1, 'hi' : 0.5},
            'b' : {'lo' : 2.0, 'hi' : 2.5},
            'i' : {'lo' : 0.5, 'hi' : 0.7},
            'x' : {'lo' : 0.7, 'hi' : 2.0}
        }
        delta = 0.05
        for char in traj_str:
            params = char_to_parameters[char]
            n_max = int(old_div((params['hi'] - params['lo']),delta))
            sequence.append(params['lo'] + delta*random.randint(1, n_max-1))

        return make_1d_traj(coordinates=sequence,
                            velocities=[1.0]*len(sequence))

    def test_analyze_continuous_time(self):
        resultA = self.analyzer.analyze_continuous_time(self.trajectory, 
                                                        self.stateA)
        assert_equal(resultA.n_frames.tolist(), [3, 1, 1, 1, 1, 1, 1])
        resultB = self.analyzer.analyze_continuous_time(self.trajectory,
                                                        self.stateB)
        assert_equal(resultB.n_frames.tolist(), [1, 1, 1, 5])
        assert_almost_equal(resultA.times.mean(), 9.0/7.0*0.1)
        assert_almost_equal(resultB.times.mean(), 8.0/4.0*0.1)

    def test_analyze_lifetime(self):
        resA = self.analyzer.analyze_lifetime(self.trajectory, self.stateA)
        resB = self.analyzer.analyze_lifetime(self.trajectory, self.stateB)
        assert_equal(resA.n_frames.tolist(), [3, 1, 2])  # A->B
        assert_equal(resB.n_frames.tolist(), [2, 1, 1, 11])  # B->A
        assert_almost_equal(resA.times.mean(), 6.0/3.0*0.1)
        assert_almost_equal(resB.times.mean(), 15.0/4.0*0.1)

    def test_analyze_transition_duration(self):
        resAB = self.analyzer.analyze_transition_duration(self.trajectory,
                                                          self.stateA,
                                                          self.stateB)
        resBA = self.analyzer.analyze_transition_duration(self.trajectory,
                                                          self.stateB,
                                                          self.stateA)
        assert_equal(resAB.n_frames.tolist(), [2, 0, 0, 1])
        assert_equal(resBA.n_frames.tolist(), [1, 0, 0, 6])
        assert_equal(resAB.times.mean(), 3.0/4.0*0.1)
        assert_equal(resBA.times.mean(), 7.0/4.0*0.1)

    def test_analyze_flux(self):
        # A: [{out: 1, in: 1}
        # B: insufficient
        # NOTE: we may want a separate trajectory for this
        flux_core_test_str = "axaxaaaxxaxbxaaxxabaa"
        # frame numbers       0    5    0    5    0
        # in (I) or out (O):   OIOIIIOOI   IIOO
        core_traj = self._make_traj(flux_core_test_str)
        flux_segs_A = self.analyzer.analyze_flux(core_traj, self.stateA)
        # flux_segs_A = self.analyzer.flux_segments[self.stateA]
        assert_equal(len(flux_segs_A['in']), 4)
        assert_equal(flux_segs_A['in'][0], core_traj[2:3])
        assert_equal(flux_segs_A['in'][1], core_traj[4:7])
        assert_equal(flux_segs_A['in'][2], core_traj[9:10])
        assert_equal(flux_segs_A['in'][3], core_traj[13:15])

        assert_equal(len(flux_segs_A['out']), 4)
        assert_equal(flux_segs_A['out'][0], core_traj[1:2])
        assert_equal(flux_segs_A['out'][1], core_traj[3:4])
        assert_equal(flux_segs_A['out'][2], core_traj[7:9])
        assert_equal(flux_segs_A['out'][3], core_traj[15:17])

    def test_analyze_flux_with_interface(self):
        flux_iface_traj_str = "aixixaiaxiixiaxaixbxbixiaaixiai"
        # frame numbers        0    5    0    5    0    5    0
        # in (I) or out (O)      OOOIIIOOOOOIOII       IIIOO
        assert_equal(len(flux_iface_traj_str), 31) # check I counted correctly
        flux_traj = self._make_traj(flux_iface_traj_str)
        self.analyzer.reset_analysis()
        flux_segs_A = self.analyzer.analyze_flux(flux_traj, self.stateA,
                                                 self.interfaceA0)
        # flux_segs_A = self.analyzer.flux_segments[self.stateA]
        assert_equal(flux_segs_A['in'][:],
                     [flux_traj[5:8], flux_traj[13:14], flux_traj[15:17],
                      flux_traj[24:27]])
        assert_equal(flux_segs_A['out'][:],
                     [flux_traj[2:5], flux_traj[8:13], flux_traj[14:15],
                      flux_traj[27:29]])

    def test_minus_flux(self):
        flux_iface_traj_str = "axxxaaaxxxa"
        flux_traj = self._make_traj(flux_iface_traj_str)
        self.analyzer.reset_analysis()
        flux_segs_A = self.analyzer.analyze_flux(flux_traj, self.stateA,
                                                 self.interfaceA0)
        assert_equal(flux_segs_A['in'][:], [flux_traj[4:7]])
        assert_equal(flux_segs_A['out'][:],
                     [flux_traj[1:4], flux_traj[7:10]])

    def test_flux(self):
        flux_iface_traj_str = "aixixaiaxiixiaxaixbxbixiaaixiai"
        flux_traj = self._make_traj(flux_iface_traj_str)
        self.analyzer.reset_analysis()
        self.analyzer.dt = 1.0
        flux = self.analyzer.flux(trajectories=[flux_traj],
                                  state=self.stateA,
                                  interface=self.interfaceA0)
        average_out = old_div((3.0 + 5.0 + 1.0 + 2.0), 4.0)
        average_in = old_div((3.0 + 1.0 + 2.0 + 3.0), 4.0)
        assert_almost_equal(flux, old_div(1.0, (average_out + average_in)))
        self.analyzer.dt = None

    @raises(RuntimeError)
    def test_flux_no_dt(self):
        analyzer = paths.TrajectoryTransitionAnalysis(self.transition)
        flux_iface_traj_str = "aixixaiaxiixiaxaixbxbixiaaixiai"
        flux_traj = self._make_traj(flux_iface_traj_str)
        flux = analyzer.flux(trajectories=[flux_traj],
                             state=self.stateA,
                             interface=self.interfaceA0)

    def test_analyze(self):
        # only test that it runs -- correctness testing in the others
        self.analyzer.analyze(self.trajectory)
        cont_frames = self.analyzer.continuous_frames
        life_frames = self.analyzer.lifetime_frames
        trans_frames = self.analyzer.transition_duration_frames
        cont_times = self.analyzer.continuous_times
        life_times = self.analyzer.lifetimes
        trans_times = self.analyzer.transition_duration

        self.analyzer.reset_analysis()
        self.analyzer.analyze([self.trajectory])
        assert_equal(cont_frames[self.stateA].tolist(),
                     self.analyzer.continuous_frames[self.stateA].tolist())
        assert_equal(life_frames[self.stateA].tolist(),
                     self.analyzer.lifetime_frames[self.stateA].tolist())
        A2B = (self.stateA, self.stateB)
        assert_equal(trans_frames[A2B].tolist(),
                     self.analyzer.transition_duration_frames[A2B].tolist())
        assert_equal(cont_times[self.stateA].mean(),
                     self.analyzer.continuous_times[self.stateA].mean())
        assert_equal(life_times[self.stateA].mean(),
                     self.analyzer.lifetimes[self.stateA].mean())
        assert_equal(trans_times[A2B].mean(),
                     self.analyzer.transition_duration[A2B].mean())

