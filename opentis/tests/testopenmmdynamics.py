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
                   'n_frames_max' : 5000,
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

    def teardown(self):
        if os.path.isfile(data_filename("openmmengine_test.nc")):
            os.remove(data_filename("openmmengine_test.nc"))

    def test_sanity(self):
        assert_equal(os.path.isfile(data_filename("openmmengine_test.nc")),
                     True)

    def test_equilibrate(self):
        assert SkipTest

    def test_snapshot_get(self):
        assert SkipTest

    def test_snapshot_set(self):
        assert SkipTest

    def test_generate_next_frame(self):
        assert SkipTest

    def test_generate(self):
        assert SkipTest
