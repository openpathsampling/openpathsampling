from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        raises, assert_almost_equal, assert_true)
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

    def test_shooting_point_analysis(self):
        raise SkipTest

    def test_committor(self):
        raise SkipTest

    def test_committor_histogram(self):
        raise SkipTest

    def test_to_pandas(self):
        raise SkipTest

