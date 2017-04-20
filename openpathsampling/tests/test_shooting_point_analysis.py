from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        raises, assert_almost_equal, assert_true, assert_in,
                        assert_raises)
from nose.plugins.skip import SkipTest
from numpy.testing import assert_array_almost_equal
from test_helpers import make_1d_traj, data_filename

import openpathsampling as paths
import openpathsampling.engines as peng
import numpy as np
import os

from openpathsampling.analysis.shooting_point_analysis import *

import logging
logging.getLogger('openpathsampling.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.netcdfplus').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.ensemble').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.engines').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.pathmover').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.sample').setLevel(logging.CRITICAL)


class testTransformedDict(object):
    def setup(self):
        self.untransformed = {(0, 1) : "a", (1, 2) : "b", (2, 3) : "c"}
        self.transformed = {0 : "a", 1 : "b", 2 : "c"}
        self.hash_function = lambda x : x[0]
        self.empty = TransformedDict(self.hash_function, {})
        self.test_dict = TransformedDict(self.hash_function,
                                         self.untransformed)

    def test_initialization(self):
        assert_equal(self.test_dict.store, self.transformed)
        assert_equal(self.test_dict.hash_representatives, 
                     {0: (0,1), 1: (1,2), 2: (2,3)})

    def test_set_get(self):
        self.empty[(5,6)] = "d"
        assert_equal(self.empty.store, {5: "d"})
        assert_equal(self.empty.hash_representatives, {5: (5,6)})
        assert_equal(self.empty[(5,6)], "d")

    def test_update(self):
        self.test_dict.update({(5,6): "d"})
        assert_equal(self.test_dict.store, 
                     {0 : "a", 1 : "b", 2 : "c", 5 : "d"})
        assert_equal(self.test_dict.hash_representatives, 
                     {0: (0,1), 1: (1,2), 2: (2,3), 5: (5,6)})

    def test_del(self):
        del self.test_dict[(0, 1)]
        assert_equal(self.test_dict.store, {1 : "b", 2 : "c"})

    def test_iter(self):
        iterated = [k for k in self.test_dict]
        for (truth, beauty) in zip(self.transformed.keys(), iterated):
            assert_equal(truth, beauty)

    def test_len(self):
        assert_equal(len(self.test_dict), 3)
        assert_equal(len(self.empty), 0)

    def test_rehash(self):
        rehashed = self.test_dict.rehash(lambda x : x[1])
        assert_equal(rehashed.store, {1: "a", 2: "b", 3: "c"})
        assert_equal(rehashed.hash_representatives,
                     {1: (0,1), 2: (1,2), 3: (2,3)})


class testSnapshotByCoordinateDict(object):
    def setup(self):
        self.empty_dict = SnapshotByCoordinateDict()
        coords_A = np.array([[0.0, 0.0]])
        coords_B = np.array([[1.0, 1.0]])
        self.key_A = coords_A.tostring()
        self.key_B = coords_B.tostring()
        self.snapA1 = peng.toy.Snapshot(coordinates=coords_A,
                                        velocities=np.array([[0.0, 0.0]]))
        self.snapA2 = peng.toy.Snapshot(coordinates=coords_A,
                                        velocities=np.array([[1.0, 1.0]]))
        self.snapB1 = peng.toy.Snapshot(coordinates=coords_B,
                                        velocities=np.array([[0.0, 0.0]]))
        self.dict1 = SnapshotByCoordinateDict({self.snapA1: "A1",
                                               self.snapB1: "B1"})

    def test_initialization(self):
        assert_equal(self.dict1.store, {self.key_A: "A1", self.key_B: "B1"})

    def test_get_set(self):
        self.dict1[self.snapA2] = "A2"
        assert_equal(self.dict1.store, {self.key_A: "A2", self.key_B: "B1"})


class testShootingPointAnalysis(object):
    def setup(self):
        # taken from the testCommittorSimulation
        import openpathsampling.engines.toy as toys
        pes = toys.LinearSlope(m=[0.0], c=[0.0]) # flat line
        topology = toys.Topology(n_spatial=1, masses=[1.0], pes=pes)
        descriptor = peng.SnapshotDescriptor.construct(
            toys.Snapshot,
            {
                'n_atoms': 1,
                'n_spatial': 1
            }
        )
        engine = peng.NoEngine(descriptor)
        self.snap0 = toys.Snapshot(coordinates=np.array([[0.0]]),
                                   velocities=np.array([[1.0]]),
                                   engine=engine)
        self.snap1 = toys.Snapshot(coordinates=np.array([[0.1]]),
                                   velocities=np.array([[1.0]]),
                                   engine=engine)
        integrator = toys.LeapfrogVerletIntegrator(0.1)
        options = {
            'integ': integrator,
            'n_frames_max': 10000,
            'n_steps_per_frame': 5
        }
        self.engine = toys.Engine(options=options, topology=topology)
        cv = paths.FunctionCV("Id", lambda snap : snap.coordinates[0][0])
        self.left = paths.CVDefinedVolume(cv, float("-inf"), -1.0)
        self.right = paths.CVDefinedVolume(cv, 1.0, float("inf"))

        randomizer = paths.NoModification()
        self.filename = data_filename("shooting_analysis.nc")
        self.storage = paths.Storage(self.filename, 
                                     mode="w")

        self.simulation = paths.CommittorSimulation(
            storage=self.storage,
            engine=self.engine,
            states=[self.left, self.right],
            randomizer=randomizer,
            initial_snapshots=[self.snap0, self.snap1]
        )
        self.simulation.output_stream = open(os.devnull, 'w')
        self.simulation.run(20)
        # set up the analysis object
        self.analyzer = ShootingPointAnalysis(self.storage.steps,
                                              [self.left, self.right])

    def teardown(self):
        import os
        if os.path.isfile(self.filename):
            os.remove(self.filename)
        paths.EngineMover.default_engine = None # set by Committor

    def test_shooting_point_analysis(self):
        assert_equal(len(self.analyzer), 2)
        assert_true(0 < self.analyzer[self.snap0][self.left] < 20)
        assert_true(0 < self.analyzer[self.snap0][self.right] < 20)
        assert_true(0 < self.analyzer[self.snap1][self.left] < 20)
        assert_true(0 < self.analyzer[self.snap1][self.right] < 20)

    def test_from_individual_runs(self):
        runs = [(self.snap0, self.left),
                (self.snap0, self.left),
                (self.snap0, self.left),
                (self.snap0, self.right),
                (self.snap1, self.left),
                (self.snap1, self.right),
                (self.snap1, self.left),
                (self.snap1, self.right)]
        analyzer = ShootingPointAnalysis.from_individual_runs(runs)
        assert_equal(analyzer[self.snap0][self.left], 3)
        assert_equal(analyzer[self.snap0][self.right], 1)
        assert_equal(analyzer[self.snap1][self.left], 2)
        assert_equal(analyzer[self.snap1][self.right], 2)

    def test_non_shooting_steps(self):
        network = paths.TPSNetwork(self.left, self.right)
        init_traj = make_1d_traj([-1.1, 0.0, 1.1])
        ensemble = network.all_ensembles[0]
        mover = paths.PathReversalMover(ensemble)
        scheme = paths.LockedMoveScheme(mover, network)
        init_conds=scheme.initial_conditions_from_trajectories([init_traj])
        assert_equal(len(init_conds), 1)
        self.storage.close()
        self.storage = paths.Storage(self.filename, "w")
        assert_equal(init_conds[ensemble].trajectory, init_traj)
        sim = paths.PathSampling(storage=self.storage,
                                 move_scheme=scheme,
                                 sample_set=init_conds)
        sim.output_stream = open(os.devnull, "w")
        sim.run(1)
        step0 = self.storage.steps[0]
        step1 = self.storage.steps[1]

        assert_equal(self.analyzer.step_key(step0), None)
        assert_equal(self.analyzer.step_key(step1), None)

        assert_equal(self.analyzer.analyze_single_step(step0), [])
        assert_equal(self.analyzer.analyze_single_step(step1), [])

    def test_committor(self):
        committor_A = self.analyzer.committor(self.left)
        committor_B = self.analyzer.committor(self.right)
        assert_true(len(committor_A) == len(committor_B) == 2)
        keys = [self.snap0, self.snap1]
        hashes0 = [self.analyzer.hash_function(k) for k in keys]
        for kA, kB in zip(committor_A.keys(), committor_B.keys()):
            hashA = self.analyzer.hash_function(kA)
            hashB = self.analyzer.hash_function(kB)
            assert_equal(hashA, hashB)
            assert_in(hashA, hashes0)
            # hash is the same; snapshot is not
        for snap in committor_A:
            assert_equal(committor_A[snap], 
                         float(self.analyzer[snap][self.left]) / 20.0)
            assert_almost_equal(committor_A[snap] + committor_B[snap], 1.0)
            assert_equal(committor_B[snap], 
                         float(self.analyzer[snap][self.right]) / 20.0)

        rehash = lambda snap : 2 * snap.xyz[0][0]
        committor_A_rehash = self.analyzer.committor(self.left, rehash)
        assert_items_equal(committor_A.values(), committor_A_rehash.values())
        for snap in committor_A.keys():
            assert_in(rehash(snap), committor_A_rehash.keys())

    def test_committor_histogram_1d(self):
        rehash = lambda snap : 2 * snap.xyz[0][0]
        input_bins = [-0.05, 0.05, 0.15, 0.25, 0.35, 0.45]
        hist, bins = self.analyzer.committor_histogram(rehash, self.left,
                                                       input_bins)
        assert_equal(len(hist), 5)
        for index in [1, 3, 4]:
            assert_true(np.isnan(hist[index]))
        for index in [0, 2]:
            assert_true(hist[index] > 0)
        assert_array_almost_equal(bins, input_bins)

    def test_committor_histogram_2d(self):
        rehash = lambda snap : (snap.xyz[0][0], 2 * snap.xyz[0][0])
        input_bins = [-0.05, 0.05, 0.15, 0.25, 0.35, 0.45]
        hist, b_x, b_y = self.analyzer.committor_histogram(rehash, self.left,
                                                           input_bins)
        assert_equal(hist.shape, (5,5))
        for i in range(5):
            for j in range(5):
                if (i,j) in [(0, 0), (1, 2)]:
                    pass
                    assert_true(hist[(i,j)] > 0)
                else:
                    assert_true(np.isnan(hist[(i,j)]))

        # this may change later to bins[0]==bins[1]==input_bins
        assert_array_almost_equal(input_bins, b_x)
        assert_array_almost_equal(input_bins, b_y)

    @raises(RuntimeError)
    def test_committor_histogram_3d(self):
        # only 1D and 2D are supported
        rehash = lambda snap : (snap.xyz[0][0], 2 * snap.xyz[0][0], 0.0)
        input_bins = [-0.05, 0.05, 0.15, 0.25, 0.35, 0.45]
        hist, b_x, b_y = self.analyzer.committor_histogram(rehash, self.left,
                                                           input_bins)

    def test_to_pandas(self):
        df1 = self.analyzer.to_pandas()
        df2 = self.analyzer.to_pandas(lambda x : x.xyz[0][0])
        assert_equal(df1.shape, (2,2))
        assert_items_equal(df1.index, range(2))
        assert_items_equal(df2.index, [0.0, 0.1])
        assert_items_equal(df1.columns, [self.left.name, self.right.name])


class testReactiveFluxAnalysis(object):
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
                                      velocities=np.array([[2.0]]),
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
                                      velocities=np.array([[3.0]]),
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
        self.analysis = ReactiveFluxAnalysis(steps=self.storage.steps,
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
                assert_almost_equal(flux_dict[snapshot], 3.0)
                assert_equal(self.analysis[snapshot]["accepted"], 2)
                assert_equal(self.analysis[snapshot]["rejected"], 1)
                assert_almost_equal(self.analysis[snapshot]["sumflux"], 6.0)
            elif snapshot.xyz[0][0] == 0.0:
                assert_almost_equal(flux_dict[snapshot], 7.0 / 4.0)
                assert_equal(self.analysis[snapshot]["accepted"], 2)
                assert_equal(self.analysis[snapshot]["rejected"], 2)
                assert_almost_equal(self.analysis[snapshot]["sumflux"], 7.0)
            elif snapshot.xyz[0][0] == 0.001:
                assert_almost_equal(flux_dict[snapshot], 22.0 / 5.0)
                assert_equal(self.analysis[snapshot]["accepted"], 4)
                assert_equal(self.analysis[snapshot]["rejected"], 1)
                assert_almost_equal(self.analysis[snapshot]["sumflux"], 22.0)

        # check total flux.
        assert_almost_equal(flux, 35.0 / 12.0)

        # check 1D histogram
        # use only two bins, so results of two snapshots are combined.
        hash1D = lambda snap: snap.xyz[0][0]
        bins1D = [-0.0015, 0.0005, 0.0015]
        hist, bins_x = self.analysis.flux_histogram(hash1D, bins1D)
        for b1, b2 in zip(bins1D, bins_x):
            assert_almost_equal(b1, b2)
        assert_almost_equal(hist[0], 13.0 / 7.0)
        assert_almost_equal(hist[1], 22.0 / 5.0)

        # check 2D histogram
        # same bins as in 1D example
        hash2D = lambda snap: (snap.xyz[0][0], snap.xyz[0][0])
        bins2D = [-0.0015, 0.0005, 0.0015]
        hist, bins_x, bins_y = self.analysis.flux_histogram(hash2D, bins2D)
        for b1, b2, b3 in zip(bins2D, bins_x, bins_y):
            assert_almost_equal(b1, b2)
            assert_almost_equal(b1, b3)
        assert_almost_equal(hist[0][0], 13.0 / 7.0)
        assert_almost_equal(hist[1][1], 22.0 / 5.0)

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
