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

    def teardown(self):
        if os.path.isfile(data_filename("storage_test.nc")):
            os.remove(data_filename("storage_test.nc"))

    def test_create_template(self):
        Storage(filename=self.filename, template=self.template_snapshot, mode='w')
        assert(os.path.isfile(data_filename("storage_test.nc")))

    def test_create_atoms(self):
        pass

    def test_stored_topology(self):
        store = Storage(filename=self.filename, template=self.template_snapshot, mode='w')
        assert(os.path.isfile(self.filename))
        store.close()

        store = Storage(filename=self.filename, mode='a')
        loaded_topology = store.template.topology

        # check if poth topologies have the same JSON string (this also tests the simplifier for topologies

        assert_equal(
            store.simplifier.topology_to_json(self.template_snapshot.topology),
            store.simplifier.topology_to_json(loaded_topology)
        )

        pass

    def test_write_load_str(self):
        store = Storage(filename=self.filename, template=self.template_snapshot, mode='w')
        assert(os.path.isfile(self.filename))

        test_str = 'test_string'
        store.init_str('test_variable')
        store.write_str('test_variable', test_str)
        store.close()

        store2 = Storage(filename=self.filename, mode='a')
        loaded_str = store2.load_str('test_variable')

        assert(loaded_str == test_str)
        pass

    def test_stored_template(self):
        store = Storage(filename=self.filename, template=self.template_snapshot, mode='w')
        assert(os.path.isfile(self.filename))
        store.close()

        store = Storage(filename=self.filename, mode='a')
        loaded_template = store.template

        compare_snapshot(loaded_template, self.template_snapshot)
        pass

    def test_load_save(self):
        store = Storage(filename=self.filename, template=self.template_snapshot, mode='w')
        assert(os.path.isfile(self.filename))

        copy = self.template_snapshot.copy()
        store.save(copy)

        store.close()

        store = Storage(filename=self.filename, mode='a')
        loaded_template = store.template

        compare_snapshot(loaded_template, self.template_snapshot)
        loaded_copy = store.load(Snapshot, 1)

        compare_snapshot(loaded_template, loaded_copy)
        pass


    def test_clone(self):
        store = Storage(filename=self.filename, template=self.template_snapshot, mode='w')
        assert(os.path.isfile(self.filename))

        copy = self.template_snapshot.copy()
        store.save(copy)

        store.save(self.traj)
        store.clone(filename=self.filename_clone, subset = self.options['solute_indices'])

        # clone the storage and reduce the number of atoms to only solute

        store2 = Storage(filename=self.filename_clone, mode='a')

        # do some tests, if this is still the same data

        compare_snapshot(
            store2.snapshot.load(0),
            store.snapshot.load(0).subset(self.options['solute_indices'])
        )

        compare_snapshot(
            store2.snapshot.load(1),
            store.snapshot.load(1).subset(self.options['solute_indices'])
        )
        store.close()
        store2.close()

        pass

    def test_clone_empty(self):
        store = Storage(filename=self.filename, template=self.template_snapshot, mode='w')
        assert(os.path.isfile(self.filename))

        copy = self.template_snapshot.copy()
        store.save(copy)

        store.save(self.traj)
        store.clone_empty(filename=self.filename_clone)

        # clone the storage and reduce the number of atoms to only solute

        store2 = Storage(filename=self.filename_clone, mode='a')

        # do some tests, if this is still the same data

        compare_snapshot(
            store2.snapshot.load(0),
            store.snapshot.load(0)
        )

        assert_equal(store2.snapshot.count(), 1)
        assert_equal(store2.trajectory.count(), 0)

        pass