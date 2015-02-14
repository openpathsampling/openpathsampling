'''
@author David W.H. Swenson
'''
import os
from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        assert_almost_equal, raises)
from nose.plugins.skip import Skip, SkipTest
from test_helpers import (true_func, data_filename,
                          assert_equal_array_array,
                          assert_not_equal_array_array)

import numpy.testing as npt

from openpathsampling.openmm_engine import *
from openpathsampling.snapshot import Snapshot
from openpathsampling.snapshot import Momentum, Configuration

import simtk.unit as u
import time

def compare_snapshot(snapshot1, snapshot2):
    npt.assert_allclose(snapshot1.box_vectors, snapshot2.box_vectors, rtol=1e-7, atol=0)
    npt.assert_allclose(snapshot1.coordinates, snapshot2.coordinates, rtol=1e-7, atol=0)
    npt.assert_allclose(snapshot1.velocities, snapshot2.velocities, rtol=1e-7, atol=0)

    assert_equal(snapshot1.potential_energy, snapshot2.potential_energy)
    assert_equal(snapshot1.kinetic_energy, snapshot2.kinetic_energy)




class testStorage(object):
    def setUp(self):
        # Use the standard Alanine to generate snapshots to store for higher testing

        self.filename = data_filename("toy_tis.nc")
        self.storage = paths.storage.Storage(filename=self.filename, mode='a')


    def teardown(self):
        pass


    def test_trajectory_iter_correct_order(self):
        traj =  self.storage.trajectory[0]

        count = 0
        for snap in traj:
            assert(traj[count] is snap)
            count += 1