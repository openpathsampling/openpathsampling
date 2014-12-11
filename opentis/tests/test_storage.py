__author__ = 'jan-hendrikprinz'

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

from opentis.openmm_engine import *
from opentis.snapshot import Snapshot
from opentis.snapshot import Momentum, Configuration
from opentis.trajectory import Trajectory

import simtk.unit as u
import time

class testStorage(object):
    def setUp(self):
        options = {'temperature' : 300.0 * u.kelvin,
                   'collision_rate' : 1.0 / u.picoseconds,
                   'timestep' : 2.0 * u.femtoseconds,
                   'nsteps_per_frame' : 10,
                   'n_frames_max' : 5,
                   'start_time' : time.time(),
                   'fn_initial_pdb' : data_filename("ala_small_traj.pdb"),
                   'platform' : 'fastest',
                   'solute_indices' : range(22),
                   'forcefield_solute' : 'amber96.xml',
                   'forcefield_solvent' : 'tip3p.xml'
                  }
        self.engine = OpenMMEngine.auto(
            filename=data_filename("storage_test.nc"),
            template=data_filename("ala_small_traj.pdb"),
            options=options
        )

        context = self.engine.simulation.context
        zero_array = np.zeros((self.engine.n_atoms, 3))
        context.setPositions(self.engine.template.coordinates)
        context.setVelocities(u.Quantity(zero_array, u.nanometers / u.picoseconds))

        self.storage = self.engine.storage

    def teardown(self):
        if os.path.isfile(data_filename("storage_test.nc")):
            os.remove(data_filename("storage_test.nc"))

    def test_sanity(self):
        assert_equal(os.path.isfile(data_filename("storage_test.nc")),
                     True)

    def test_template(self):
        template = self.storage.template
        snap0 = self.storage.snapshot.load(0)

        assert_equal(template, snap0)

    def test_trajectory(self):
        self.engine.current_snapshot = self.storage.template

        traj1 = Trajectory([ self.engine.generate_next_frame() for i in range(20) ])
        traj2 = Trajectory([ self.engine.generate_next_frame() for i in range(10) ])
        traj3 = traj1[0:10] + traj2[0:5]

        self.storage.save(traj1)
        self.storage.save(traj2)
        self.storage.save(traj3)

        traj4 = self.storage.load('trajectory', 0)

        # Check if caching works
        assert_equal(traj1, traj4)

        self.storage.trajectory.flush_cache()
        traj5 = self.storage.load(Trajectory, 0)

        # after flush traj should not be equal
        assert_equal(traj1 is traj5, False)

        # but all of their snapshots since they are still cached
        assert_equal(len(traj1),len(traj5))
        for no in range(len(traj1)):
            assert_equal(traj1[no], traj5[no])
