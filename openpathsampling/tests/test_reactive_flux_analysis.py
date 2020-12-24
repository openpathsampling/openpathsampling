from __future__ import division
from __future__ import absolute_import
from builtins import zip
from builtins import range
from past.utils import old_div
from builtins import object
from nose.tools import (assert_equal, assert_not_equal, raises,
                        assert_almost_equal, assert_true, assert_in,
                        assert_raises)
from nose.plugins.skip import SkipTest
from numpy.testing import assert_array_almost_equal
from .test_helpers import (make_1d_traj, data_filename, assert_items_equal,
                           assert_same_items)

import openpathsampling as paths
import openpathsampling.engines as peng
import numpy as np
import os

#from openpathsampling.analysis.shooting_point_analysis import *

import logging
logging.getLogger('openpathsampling.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.netcdfplus').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.ensemble').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.engines').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.pathmover').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.sample').setLevel(logging.CRITICAL)


class TestReactiveFluxAnalysis(object):
    def setup(self):
        import openpathsampling.engines.toy as toys
        # PES is one-dimensional linear slope (y(x) = -x)
        pes = toys.LinearSlope(m=[-1.0], c=[0.0])
        # one particle with mass 1.0
        topology = toys.Topology(n_spatial=1, masses=[1.0], pes=pes)
        integrator = toys.LeapfrogVerletIntegrator(0.02)
        options = {
            'integ' : integrator,
            'n_frames_max' : 1000,
            'n_steps_per_frame' : 5
        }
        self.engine = toys.Engine(options=options, topology=topology)
        # test uses snapshots with different velocities and slightly
        # different positions (0, -0.001, 0.001).
        # 0, 3, 6: direction ok, velocity too low => falls back to dividing
        #          surface
        # 1, 4, 7: wrong direction => backward shot towards B
        # 2, 5, 8, 9, 10, 11: direction ok, velocity high enough => successful
        #                     new trajectory
        self.initial_snapshots = [toys.Snapshot(
                                      coordinates=np.array([[-0.001]]),
                                      velocities=np.array([[1.0]]),
                                      engine=self.engine),
                                  toys.Snapshot(
                                      coordinates=np.array([[-0.001]]),
                                      velocities=np.array([[-1.0]]),
                                      engine=self.engine),
                                  toys.Snapshot(
                                      coordinates=np.array([[-0.001]]),
                                      velocities=np.array([[4.0]]),
                                      engine=self.engine),
                                  toys.Snapshot(
                                      coordinates=np.array([[0.0]]),
                                      velocities=np.array([[1.0]]),
                                      engine=self.engine),
                                  toys.Snapshot(
                                      coordinates=np.array([[0.0]]),
                                      velocities=np.array([[-1.0]]),
                                      engine=self.engine),
                                  toys.Snapshot(
                                      coordinates=np.array([[0.0]]),
                                      velocities=np.array([[2.0]]),
                                      engine=self.engine),
                                  toys.Snapshot(
                                      coordinates=np.array([[0.001]]),
                                      velocities=np.array([[1.0]]),
                                      engine=self.engine),
                                  toys.Snapshot(
                                      coordinates=np.array([[0.001]]),
                                      velocities=np.array([[-1.0]]),
                                      engine=self.engine),
                                  toys.Snapshot(
                                      coordinates=np.array([[0.001]]),
                                      velocities=np.array([[8.0]]),
                                      engine=self.engine),
                                  toys.Snapshot(
                                      coordinates=np.array([[0.0]]),
                                      velocities=np.array([[5.0]]),
                                      engine=self.engine),
                                  toys.Snapshot(
                                      coordinates=np.array([[0.001]]),
                                      velocities=np.array([[2.0]]),
                                      engine=self.engine),
                                  toys.Snapshot(
                                      coordinates=np.array([[0.001]]),
                                      velocities=np.array([[9.0]]),
                                      engine=self.engine)]
        # reaction coordinate is just x coordinate
        rc = paths.FunctionCV("Id", lambda snap : snap.coordinates[0][0])
        # state A: [-inf, -1]
        self.state_A = paths.CVDefinedVolume(rc, float("-inf"), -1.0)
        # area between A and dividing surface: [-1, 0]
        self.towards_A = paths.CVDefinedVolume(rc, -1.0, 0.0)
        # state B: [1, inf]
        self.state_B = paths.CVDefinedVolume(rc, 1.0, float("inf"))
        # define state labels
        self.state_labels = {
            "A" : self.state_A,
            "B" : self.state_B,
            "ToA": self.towards_A,
            "None" :~(self.state_A | self.state_B | self.towards_A)}

        # velocities are not randomized
        randomizer = paths.NoModification()

        self.filename = data_filename("rf_test.nc")
        self.storage = paths.Storage(self.filename, mode="w")
        self.storage.save(self.initial_snapshots)

        self.simulation = paths.ReactiveFluxSimulation(
                              storage=self.storage,
                              engine=self.engine,
                              states=[self.state_A, self.state_B],
                              randomizer=randomizer,
                              initial_snapshots=self.initial_snapshots,
                              rc=rc)
        self.simulation.output_stream = open(os.devnull, 'w')
        self.simulation.run(n_per_snapshot=1)
        self.analysis = None

    def teardown(self):
        if os.path.isfile(self.filename):
            os.remove(self.filename)
        paths.EngineMover.default_engine = None

    def gradient(self, snapshot):
        assert_equal(len(snapshot.xyz), 1)
        assert_equal(len(snapshot.xyz[0]), 1)
        return np.array([[1.0]])

    def test_analysis(self):
        self.storage = paths.Storage(self.filename, mode="r")
        self.analysis = paths.ReactiveFluxAnalysis(steps=self.storage.steps,
                                                   gradient=self.gradient)

        # check wrong analyze_single_step() argument
        empty_list = self.analysis.analyze_single_step(3.1415)
        assert_equal(empty_list, [])

        # dictionary with three entries returned (mind: trajectories start from
        # x = -0.001, 0.0, 0.001).
        assert_equal(len(self.analysis), 3)

        # get total and per-snapshot flux.
        flux, flux_dict = self.analysis.flux()

        # analyze counters
        for snapshot in flux_dict.keys():
            if snapshot.xyz[0][0] == -0.001:
                assert_almost_equal(flux_dict[snapshot], 4.0 / 3.0)
                assert_equal(self.analysis[snapshot]["accepted"], 1)
                assert_equal(self.analysis[snapshot]["rejected"], 2)
                assert_almost_equal(self.analysis[snapshot]["sumflux"], 4.0)
            elif snapshot.xyz[0][0] == 0.0:
                assert_almost_equal(flux_dict[snapshot], 7.0 / 4.0)
                assert_equal(self.analysis[snapshot]["accepted"], 2)
                assert_equal(self.analysis[snapshot]["rejected"], 2)
                assert_almost_equal(self.analysis[snapshot]["sumflux"], 7.0)
            elif snapshot.xyz[0][0] == 0.001:
                assert_almost_equal(flux_dict[snapshot], 19.0 / 5.0)
                assert_equal(self.analysis[snapshot]["accepted"], 3)
                assert_equal(self.analysis[snapshot]["rejected"], 2)
                assert_almost_equal(self.analysis[snapshot]["sumflux"], 19.0)

        # check total flux.
        assert_almost_equal(flux, 30.0 / 12.0)

        # check 1D histogram
        # use only two bins, so results of two snapshots are combined.
        hash1D = lambda snap: snap.xyz[0][0]
        bins1D = [-0.0015, 0.0005, 0.0015]
        hist, bins_x = self.analysis.flux_histogram(hash1D, bins1D)
        for b1, b2 in zip(bins1D, bins_x):
            assert_almost_equal(b1, b2)
        assert_almost_equal(hist[0], 11.0 / 7.0)
        assert_almost_equal(hist[1], 19.0 / 5.0)

        # check 2D histogram
        # same bins as in 1D example
        hash2D = lambda snap: (snap.xyz[0][0], snap.xyz[0][0])
        bins2D = [-0.0015, 0.0005, 0.0015]
        hist, bins_x, bins_y = self.analysis.flux_histogram(hash2D, bins2D)
        for b1, b2, b3 in zip(bins2D, bins_x, bins_y):
            assert_almost_equal(b1, b2)
            assert_almost_equal(b1, b3)
        assert_almost_equal(hist[0][0], 11.0 / 7.0)
        assert_almost_equal(hist[1][1], 19.0 / 5.0)

        # check failure of 3D histogram
        hash3D = lambda snap: (snap.xyz[0][0], snap.xyz[0][0], snap.xyz[0][0])
        bins3D = [-0.0015, 0.0005, 0.0015]
        assert_raises(RuntimeError,
                      self.analysis.flux_histogram,
                      hash3D,
                      bins3D)

        assert_raises(NotImplementedError, self.analysis.committor)
        assert_raises(NotImplementedError,
                      self.analysis.committor_histogram,
                      hash1D,
                      bins1D)
