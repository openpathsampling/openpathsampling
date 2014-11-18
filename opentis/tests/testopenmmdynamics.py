'''
@author David W.H. Swenson
'''
import os
from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        assert_almost_equal, raises)
from nose.plugins.skip import Skip, SkipTest
from test_helpers import true_func, data_filename

from opentis.openmm_engine import *
from opentis.snapshot import Snapshot, Momentum, Configuration

from simtk.unit import femtoseconds, picoseconds, nanometers, kelvin, dalton
from simtk.unit import Quantity
import time

class testOpenMMEngine(object):
    def setUp(self):
        options = {'temperature' : 300.0 * kelvin,
                   'collision_rate' : 1.0 / picoseconds,
                   'timestep' : 2.0 * femtoseconds,
                   'nsteps_per_frame' : 10,
                   'n_frames_max' : 5,
                   'start_time' : time.time(),
                   'fn_initial_pdb' : data_filename("ala_small_traj.pdb"),
                   'platform' : 'fastest',
                   'solute_indices' : range(22), 
                   'forcefield_solute' : 'amber96.xml',
                   'forcefield_solvent' : 'tip3p.xml'
                  }
        self.engine = OpenMMEngine(
            filename=data_filename("openmmengine_test.nc"), 
            topology_file=data_filename("ala_small_traj.pdb"), 
            opts=options, 
            mode='create'
        )
        self.engine.simulation.context.setPositions(self.engine.pdb.positions)

    def teardown(self):
        if os.path.isfile(data_filename("openmmengine_test.nc")):
            os.remove(data_filename("openmmengine_test.nc"))

    def test_sanity(self):
        assert_equal(os.path.isfile(data_filename("openmmengine_test.nc")),
                     True)

    def test_equilibrate(self):
        snap = self.engine.current_snapshot
        self.engine.equilibrate(5)
        # TODO: assertion
        assert SkipTest

    def test_snapshot_get(self):
        # TODO assertion
        assert SkipTest

    def test_snapshot_set(self):
        pdb_pos = (self.engine.pdb.positions / nanometers)
        testvel = []
        testpos = []
        for i in range(len(pdb_pos)):
            testpos.append(np.array(pdb_pos[i]) + np.array([1.0, 1.0, 1.0]))
            testvel.append([0.1*i, 0.1*i, 0.1*i])

        self.engine.current_snapshot = Snapshot(
            coordinates=testpos,
            velocities=testvel
        )
        # TODO assertion
        assert SkipTest

    def test_generate_next_frame(self):
        # test for success, not accuracy
        # assert that the snapshot has changed
        assert SkipTest

    def test_generate(self):
        # TODO; set initial conditions
        self.engine.initialized = True
        traj = self.engine.generate(self.engine.current_snapshot, [true_func])
        assert_equal(len(traj), self.engine.n_frames_max)
