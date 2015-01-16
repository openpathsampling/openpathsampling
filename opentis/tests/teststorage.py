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

import simtk.unit as u
import time



class testStorage(object):
    def setUp(self):
        # Use the standard Alanine to generate snapshots to store for higher testing

        self.options = {'temperature' : 300.0 * u.kelvin,
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

        # create a template snapshot
        self.template_snapshot = ops.snapshot_from_pdb(data_filename("ala_small_traj.pdb"))

        # and an openmm engine
        self.engine = ops.OpenMMEngine(options=self.options, template=self.template_snapshot)
        self.engine.initialized = True

        # run a small trajectory of a few steps that can be used to save, etc...
        self.traj = self.engine.generate(self.template_snapshot, running=[ops.LengthEnsemble(10).can_append])

    def teardown(self):
        if os.path.isfile(data_filename("storage_test.nc")):
            os.remove(data_filename("storage_test.nc"))

    def test_create_template(self):
        pass

    def test_create_atoms(self):
        pass

    def test_stored_topology(self):
        pass

    def test_stored_template(self):
        pass

    def test_write_str(self):
        pass

    def test_init_str(self):
        pass

    def test_save(self):
        pass

    def test_load(self):
        pass

    def test_clone(self):
        pass

    def test_clone_empty(self):
        pass

    def test_clone_storage(self):
        pass
