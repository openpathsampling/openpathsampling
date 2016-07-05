'''
@author David W.H. Swenson
'''
import os

import mdtraj as md
from nose.tools import (assert_equal)

import openpathsampling as paths

import openpathsampling.engines.openmm as peng
import openpathsampling.engines.toy as toys

from openpathsampling.netcdfplus import ObjectJSON
from openpathsampling.storage import Storage
from test_helpers import (data_filename,
                          compare_snapshot
                          )

import numpy as np
from nose.plugins.skip import SkipTest


class testStorage(object):
    def setUp(self):
        self.mdtraj = md.load(data_filename("ala_small_traj.pdb"))
        self.traj = peng.trajectory_from_mdtraj(self.mdtraj, simple_topology=True)

        self.filename = data_filename("storage_test.nc")
        self.filename_clone = data_filename("storage_test_clone.nc")

        self.simplifier = ObjectJSON()
        self.template_snapshot = self.traj[0]
        self.solute_indices = range(22)

        self.toy_topology = toys.Topology(
            n_spatial=2,
            masses=[1.0, 1.0],
            pes=None
        )

        self.engine = toys.Engine({}, self.toy_topology)

        self.toy_template = toys.Snapshot(
            coordinates=np.array([[-0.5, -0.5]]),
            velocities=np.array([[0.0,0.0]]),
            engine=self.engine
        )

    def teardown(self):
        if os.path.isfile(self.filename):
            os.remove(self.filename)

        if os.path.isfile(self.filename_clone):
            os.remove(self.filename_clone)

    def test_create_storage(self):
        store = Storage(filename=self.filename, mode='w')
        assert(os.path.isfile(data_filename("storage_test.nc")))
        store.close()

    def test_stored_topology(self):
        raise SkipTest
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
        raise SkipTest
        store = Storage(filename=self.filename, template=self.template_snapshot, mode='w')
        assert(os.path.isfile(self.filename))
        store.close()

        store = Storage(filename=self.filename, mode='a')
        loaded_template = store.template

        compare_snapshot(loaded_template, self.template_snapshot, True)

        store.close()

    def test_load_save(self):
        store = Storage(filename=self.filename, mode='w', use_uuid=False)
        assert(os.path.isfile(self.filename))

        store.save(self.template_snapshot)
        store.close()

        store = Storage(filename=self.filename, mode='a')
        loaded_template = store.snapshots[0]
        loaded_r = store.snapshots[1]

        compare_snapshot(loaded_template, self.template_snapshot, True)
        compare_snapshot(loaded_template.reversed, self.template_snapshot.reversed, True)
        compare_snapshot(loaded_r, self.template_snapshot.reversed)

        store.close()

    def test_load_save_uuid(self):
        store = Storage(filename=self.filename, mode='w')
        assert(os.path.isfile(self.filename))

        store.save(self.template_snapshot)
        store.close()

        store = Storage(filename=self.filename, mode='a')
        loaded_template = store.snapshots[self.template_snapshot.__uuid__]
        loaded_r = store.snapshots[self.template_snapshot.reversed.__uuid__]

        compare_snapshot(loaded_template, self.template_snapshot, True)
        compare_snapshot(loaded_template.reversed, self.template_snapshot.reversed, True)
        compare_snapshot(loaded_r, self.template_snapshot.reversed)

        store.close()

    def test_load_save_toy(self):
        store = Storage(filename=self.filename, mode='w', use_uuid=False)
        assert(os.path.isfile(self.filename))

        store.save(self.toy_template)

        store.close()

        store = Storage(filename=self.filename, mode='a')
        loaded_template = store.snapshots[0]
        loaded_r = store.snapshots[1]

        compare_snapshot(loaded_template, self.toy_template, True)
        compare_snapshot(loaded_template.reversed, self.toy_template.reversed, True)
        compare_snapshot(loaded_r, self.toy_template.reversed)

        store.close()

    def test_clone(self):
        raise SkipTest
        store = Storage(filename=self.filename, template=self.template_snapshot, mode='w')
        assert(os.path.isfile(self.filename))

        copy = self.template_snapshot.copy()
        store.save(copy)

        store.save(self.traj)

        store.clone(filename=self.filename_clone)

        # clone the storage and reduce the number of atoms to only solute

        store2 = Storage(filename=self.filename_clone, mode='a')

        # do some tests, if this is still the same data

        compare_snapshot(
            store2.snapshots.load(0),
            store.snapshots.load(0),
            True
        )

        compare_snapshot(
            store2.snapshots.load(1),
            store.snapshots.load(1),
            True
        )
        store.close()
        store2.close()

    def test_clone_empty(self):
        raise SkipTest

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
            store.snapshots.load(0),
            True
        )

        # check if the reversed copy also works
        compare_snapshot(
            store2.snapshots.load(1),
            store.snapshots.load(1),
            True
        )

        assert_equal(len(store2.snapshots), 2)
        assert_equal(len(store2.trajectories), 0)

        store.close()
        store2.close()

    def test_reverse_bug(self):
        store = Storage(filename=self.filename,
                        mode='w', use_uuid=False)
        assert(os.path.isfile(self.filename))

        store.snapshots.save(self.template_snapshot)
        rev = self.template_snapshot.reversed

        # save the reversed one
        store.snapshots.save(rev)

        # check that the reversed one has index 1 and not 3!
        assert(store.idx(rev) == 1)

        # and we have exactly one snapshot
        assert(len(store.snapshots) == 2)
        assert(len(store.dimensions['snapshots']) == 1)
        store.close()

    def test_version(self):
        store = Storage(filename=self.filename, template=self.template_snapshot, mode='w')
        assert(os.path.isfile(self.filename))

        assert(store.storage_version == paths.version.version)
