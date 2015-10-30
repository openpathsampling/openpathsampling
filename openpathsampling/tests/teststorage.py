'''
@author David W.H. Swenson
'''
import os
from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        assert_almost_equal, raises)
from nose.plugins.skip import Skip, SkipTest
from test_helpers import (true_func,
                          data_filename,
                          assert_equal_array_array,
                          assert_not_equal_array_array,
                          compare_snapshot
                          )

from openpathsampling.openmm_engine import *
from openpathsampling.snapshot import Snapshot
from openpathsampling.storage.object_json import ObjectJSON

import mdtraj as md


class testStorage(object):
    def setUp(self):
        self.mdtraj = md.load(data_filename("ala_small_traj.pdb"))
        self.traj = paths.tools.trajectory_from_mdtraj(self.mdtraj)

        self.filename = data_filename("storage_test.nc")
        self.filename_clone = data_filename("storage_test_clone.nc")

        self.simplifier = ObjectJSON()
        self.template_snapshot = self.traj[0]
        self.solute_indices = range(22)

    def teardown(self):
        if os.path.isfile(self.filename):
            os.remove(self.filename)

        if os.path.isfile(self.filename_clone):
            os.remove(self.filename_clone)

    def test_create_template(self):
        store = Storage(filename=self.filename, template=self.template_snapshot, mode='w')
        assert(os.path.isfile(data_filename("storage_test.nc")))
        store.close()

    def test_stored_topology(self):
        store = Storage(filename=self.filename, template=self.template_snapshot, mode='w')
        assert(os.path.isfile(self.filename))
        store.close()

        store = Storage(filename=self.filename, mode='a')
        loaded_topology = store.template.topology

        # check if poth topologies have the same JSON string (this also tests the simplifier for topologies

        assert_equal(
            self.simplifier.to_json(self.template_snapshot.topology),
            self.simplifier.to_json(loaded_topology)
        )

        store.close()

    def test_stored_template(self):
        store = Storage(filename=self.filename, template=self.template_snapshot, mode='w')
        assert(os.path.isfile(self.filename))
        store.close()

        store = Storage(filename=self.filename, mode='a')
        loaded_template = store.template

        compare_snapshot(loaded_template, self.template_snapshot)

        store.close()

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

        store.close()

    def test_clone(self):
        store = Storage(filename=self.filename, template=self.template_snapshot, mode='w')
        assert(os.path.isfile(self.filename))

        copy = self.template_snapshot.copy()
        store.save(copy)

        store.save(self.traj)

        store.clone(filename=self.filename_clone, subset = self.solute_indices)

        # clone the storage and reduce the number of atoms to only solute

        store2 = Storage(filename=self.filename_clone, mode='a')

        # do some tests, if this is still the same data

        compare_snapshot(
            store2.snapshots.load(0),
            store.snapshots.load(0).subset(self.solute_indices)
        )

        compare_snapshot(
            store2.snapshots.load(1),
            store.snapshots.load(1).subset(self.solute_indices)
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
            store2.snapshots.load(0),
            store.snapshots.load(0)
        )

        # check if the reversed copy also works
        compare_snapshot(
            store2.snapshots.load(1),
            store.snapshots.load(1)
        )

        assert_equal(len(store2.snapshots), 2)
        assert_equal(len(store2.trajectories), 0)

        store.close()
        store2.close()

        pass
