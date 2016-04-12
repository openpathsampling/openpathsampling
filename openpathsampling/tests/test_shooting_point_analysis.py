from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        raises, assert_almost_equal, assert_true)
from nose.plugins.skip import SkipTest
from numpy.testing import assert_array_almost_equal
from test_helpers import make_1d_traj

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
        pass

    def test_shooting_point_analysis(self):
        raise SkipTest
