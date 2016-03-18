from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        raises, assert_almost_equal, assert_true)
from nose.plugins.skip import SkipTest
from numpy.testing import assert_array_almost_equal
from test_helpers import make_1d_traj

import openpathsampling as paths
import openpathsampling.engines as peng
import numpy as np

from openpathsampling.snapshot_modifier import *

import logging
logging.getLogger('openpathsampling.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.netcdfplus').setLevel(logging.CRITICAL)

class testNoModification(object):
    def setup(self):
        self.modifier = NoModification()
        self.snapshot_1D = peng.toy.Snapshot(
            coordinates=np.array([0.0, 1.0, 2.0, 3.0]), 
            velocities=np.array([0.5, 1.5, 2.5, 3.5])
        )
        self.snapshot_3D = peng.openmm.MDSnapshot(
            coordinates=np.array([[0.0, 0.1, 0.2], 
                                  [1.0, 1.1, 1.2],
                                  [2.0, 2.1, 2.2],
                                  [3.0, 3.1, 3.2]]),
            velocities=np.array([[0.5, 0.6, 0.7],
                                 [1.5, 1.6, 1.7],
                                 [2.5, 2.6, 2.7],
                                 [3.5, 3.6, 3.7]])
        )

    # test methods from the abstract base class along with NoModification
    def test_extract_subset(self):
        mod = NoModification(subset_mask=[1,2])
        sub_1Dx = mod.extract_subset(self.snapshot_1D.coordinates)
        assert_array_almost_equal(sub_1Dx, np.array([1.0, 2.0]))
        sub_1Dv = mod.extract_subset(self.snapshot_1D.velocities)
        assert_array_almost_equal(sub_1Dv, np.array([1.5, 2.5]))

        sub_3Dx = mod.extract_subset(self.snapshot_3D.coordinates)
        assert_array_almost_equal(sub_3Dx, np.array([[1.0, 1.1, 1.2],
                                                     [2.0, 2.1, 2.2]]))
        sub_3Dv = mod.extract_subset(self.snapshot_3D.velocities)
        assert_array_almost_equal(sub_3Dv, np.array([[1.5, 1.6, 1.7],
                                                     [2.5, 2.6, 2.7]]))

    def test_apply_to_subset(self):
        mod = NoModification(subset_mask=[1,2])
        copy_1Dx = self.snapshot_1D.coordinates.copy()
        new_1Dx = mod.apply_to_subset(copy_1Dx, np.array([-1.0, -2.0]))
        assert_array_almost_equal(new_1Dx, np.array([0.0, -1.0, -2.0, 3.0]))
        # and check that memory points to the right things; orig unchanged
        assert_true(copy_1Dx is new_1Dx)
        assert_array_almost_equal(self.snapshot_1D.coordinates,
                                  np.array([0.0, 1.0, 2.0, 3.0]))

        copy_3Dx = self.snapshot_3D.coordinates.copy()
        new_3Dx = mod.apply_to_subset(copy_3Dx, 
                                      np.array([[-1.0, -1.1, -1.2], 
                                                [-2.0, -2.1, -2.2]]))
        assert_array_almost_equal(new_3Dx, np.array([[0.0, 0.1, 0.2], 
                                                     [-1.0, -1.1, -1.2],
                                                     [-2.0, -2.1, -2.2],
                                                     [3.0, 3.1, 3.2]]))
        # and check that memory points to the right things; orig unchanged
        assert_true(copy_3Dx is new_3Dx)
        assert_array_almost_equal(self.snapshot_3D.coordinates,
                                  np.array([[0.0, 0.1, 0.2], 
                                            [1.0, 1.1, 1.2],
                                            [2.0, 2.1, 2.2],
                                            [3.0, 3.1, 3.2]]))

    def test_call(self):
        new_1D = self.modifier(self.snapshot_1D)
        assert_array_almost_equal(self.snapshot_1D.coordinates,
                                  new_1D.coordinates)
        assert_array_almost_equal(self.snapshot_1D.velocities,
                                  new_1D.velocities)
        new_3D = self.modifier(self.snapshot_3D)
        assert_array_almost_equal(self.snapshot_3D.coordinates,
                                  new_3D.coordinates)
        assert_array_almost_equal(self.snapshot_3D.velocities,
                                  new_3D.velocities)
        # TODO: this is awaiting changes in #445
        # assert_true(self.snapshot_1D.coordinates is not new_1D.coordinates)
        # assert_true(self.snapshot_1D.velocities is not new_1D.velocities)
        # assert_true(self.snapshot_3D.coordinates is not new_3D.coordinates)
        # assert_true(self.snapshot_3D.velocities is not new_3D.velocities)


class testRandomizeVelocities(object):
    def setup(self):
        pass

    def test_call(self):
        raise SkipTest

    def test_subset_call(self):
        raise SkipTest
