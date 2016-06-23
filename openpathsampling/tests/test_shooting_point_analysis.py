from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        raises, assert_almost_equal, assert_true, assert_in)
from nose.plugins.skip import SkipTest
from numpy.testing import assert_array_almost_equal
from test_helpers import make_1d_traj, data_filename

import openpathsampling as paths
import openpathsampling.engines as peng
import numpy as np

from openpathsampling.analysis.shooting_point_analysis import *

import logging
logging.getLogger('openpathsampling.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.netcdfplus').setLevel(logging.CRITICAL)

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
        self.snap0 = toys.Snapshot(coordinates=np.array([[0.0]]),
                                   velocities=np.array([[1.0]]),
                                   topology=topology)
        self.snap1 = toys.Snapshot(coordinates=np.array([[0.1]]),
                                   velocities=np.array([[1.0]]),
                                   topology=topology)
        integrator = toys.LeapfrogVerletIntegrator(0.1)
        options = {
            'integ': integrator,
            'n_frames_max': 10000,
            'nsteps_per_frame': 5
        }
        self.engine = toys.Engine(options=options, template=self.snap0)
        cv = paths.CV_Function("Id", lambda snap : snap.coordinates[0][0])
        self.left = paths.CVRangeVolume(cv, float("-inf"), -1.0)
        self.right = paths.CVRangeVolume(cv, 1.0, float("inf"))

        randomizer = paths.NoModification()
        self.filename = data_filename("shooting_analysis.nc")
        self.storage = paths.Storage(self.filename, 
                                     mode="w", 
                                     template=self.snap0)

        self.simulation = paths.CommittorSimulation(
            storage=self.storage,
            engine=self.engine,
            states=[self.left, self.right],
            randomizer=randomizer,
            initial_snapshots=[self.snap0, self.snap1]
        )
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
        assert_equal(bins, input_bins)

    def test_committor_histogram_2d(self):
        rehash = lambda snap : (snap.xyz[0][0], 2 * snap.xyz[0][0])
        input_bins = [-0.05, 0.05, 0.15, 0.25, 0.35, 0.45]
        hist, bins = self.analyzer.committor_histogram(rehash, self.left,
                                                       input_bins)
        assert_equal(hist.shape, (5,5))
        for i in range(5):
            for j in range(5):
                if (i,j) in [(0, 0), (1, 2)]:
                    assert_true(hist[(i,j)] > 0)
                else:
                    assert_true(np.isnan(hist[(i,j)]))

        # this may change later to bins[0]==bins[1]==input_bins
        assert_equal(bins, input_bins) 

    def test_to_pandas(self):
        df1 = self.analyzer.to_pandas()
        df2 = self.analyzer.to_pandas(lambda x : x.xyz[0][0])
        assert_equal(df1.shape, (2,2))
        assert_items_equal(df1.index, range(2))
        assert_items_equal(df2.index, [0.0, 0.1])
        assert_items_equal(df1.columns, [self.left.name, self.right.name])

