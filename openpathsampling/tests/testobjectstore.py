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

import copy

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
        self.template_snapshot = paths.snapshot_from_pdb(data_filename("ala_small_traj.pdb"))

        # and an openmm engine
        self.engine = paths.OpenMMEngine(options=self.options, template=self.template_snapshot)
        self.engine.initialized = True

        # run a small trajectory of a few steps that can be used to save, etc...
        self.traj = self.engine.generate(self.template_snapshot, running=[paths.LengthEnsemble(2).can_append])

        self.filename = data_filename("storage_test.nc")
        self.filename_clone = data_filename("storage_test_clone.nc")
        self.filename_demo = data_filename("toy_tis.nc")

    def teardown(self):
        if os.path.isfile(self.filename):
            os.remove(self.filename)

#        if os.path.isfile(self.filename_clone):
#            os.remove(self.filename_clone)

    def test_store_orderparameters(self):
        storage = Storage(filename=self.filename_demo, mode='a')
        storage.clone(filename=self.filename_clone)
        assert(os.path.isfile(self.filename_clone))

        storage2 = Storage(filename=self.filename_clone, mode='a')
        for store_name in storage.list_stores():
            print '%s' % store_name
            store = getattr(storage, store_name)

            if hasattr(store.content_class, 'creatable'):
                for obj in store:
                    storage2.save(obj)

        storage.close()
        storage2.close()