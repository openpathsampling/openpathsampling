"""
@author David W.H. Swenson
@author Jan-Hendrik Prinz
"""
from __future__ import absolute_import
from builtins import zip
from builtins import range
from builtins import object
import os

import pytest

from nose.tools import (assert_equal)

import openpathsampling as paths

import openpathsampling.engines.openmm as peng
import openpathsampling.engines.toy as toys

from openpathsampling.netcdfplus import ObjectJSON
from openpathsampling.storage import Storage
from .test_helpers import (data_filename, md, compare_snapshot)

import numpy as np
from nose.plugins.skip import SkipTest


class TestStorage(object):
    def setup(self):
        if not md:
            raise SkipTest("mdtraj not installed")
        self.mdtraj = md.load(data_filename("ala_small_traj.pdb"))
        _ = pytest.importorskip('simtk.unit')
        self.traj = peng.trajectory_from_mdtraj(
            self.mdtraj, simple_topology=True)

        self.filename = data_filename("storage_test.nc")
        self.filename_clone = data_filename("storage_test_clone.nc")

        self.simplifier = ObjectJSON()
        self.template_snapshot = self.traj[0]
        self.solute_indices = list(range(22))

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
        store = Storage(
            filename=self.filename,
            mode='w')
        assert(os.path.isfile(self.filename))
        store.close()

        store = Storage(filename=self.filename, mode='a')
        loaded_topology = store.template.topology
        # check if path topologies have the same JSON string
        # this also tests the simplifier for topologies

        assert_equal(
            self.simplifier.to_json(self.template_snapshot.topology),
            self.simplifier.to_json(loaded_topology)
        )

        store.close()

    def test_safemode(self):
        fname = data_filename("cv_storage_safemode_test.nc")
        if os.path.isfile(fname):
            os.remove(fname)

        cv = paths.CoordinateFunctionCV('cv', lambda x: x)

        traj = paths.Trajectory(list(self.traj))
        template = traj[0]

        storage_w = paths.Storage(fname, "w")
        storage_w.snapshots.save(template)
        storage_w.cvs.save(cv)
        storage_w.close()

        storage_r = paths.Storage(fname, 'r')
        # default safemode = False
        assert(storage_r.simplifier.safemode is False)
        cv_r = storage_r.cvs[0]
        assert(cv_r == cv)
        assert(cv.cv_callable is not None)
        storage_r.close()

        storage_r = paths.Storage(fname, 'r')
        storage_r.simplifier.safemode = True
        cv_r = storage_r.cvs[0]
        assert(cv_r == cv)
        assert(cv_r.cv_callable is None)
        storage_r.close()

    def test_store_snapshots(self):
        fname = data_filename("cv_storage_test.nc")
        if os.path.isfile(fname):
            os.remove(fname)

        traj = paths.Trajectory(list(self.traj))
        template = traj[0]

        for use_cache in (False, True):
            # print '=========================================================='
            # print 'UUID', use_uuid, 'CACHE', use_cache
            # print '=========================================================='

            storage_w = paths.Storage(fname, "w")
            storage_w.snapshots.save(template)

            # let's mess up the order in which we save and include
            # reversed ones as well

            assert(len(storage_w.snapshots) == 2)
            assert(len(storage_w.trajectories) == 0)
            assert(len(storage_w.stores['snapshot0']) == 2)
            storage_w.snapshots.save(traj[8].reversed)
            assert(len(storage_w.snapshots) == 4)
            assert(len(storage_w.trajectories) == 0)
            assert(len(storage_w.stores['snapshot0']) == 4)
            # this will store traj[6:] under pos IDX #0
            storage_w.trajectories.save(traj[6:])
            assert(len(storage_w.snapshots) == 10)
            assert(len(storage_w.trajectories) == 1)
            assert(len(storage_w.stores['snapshot0']) == 10)

            traj_rev = traj.reversed

            # this will store traj_rev under pos IDX #1
            storage_w.trajectories.mention(traj_rev)
            assert(len(storage_w.snapshots) == 20)
            assert(len(storage_w.trajectories) == 2)
            assert(len(storage_w.stores['snapshot0']) == 10)

            # this will not do anything since traj is already saved
            storage_w.trajectories.save(traj_rev)
            assert(len(storage_w.snapshots) == 20)
            assert(len(storage_w.trajectories) == 2)
            assert(len(storage_w.stores['snapshot0']) == 10)

            # this will store traj under pos IDX #2
            storage_w.trajectories.save(traj)
            assert(len(storage_w.snapshots) == 20)
            assert(len(storage_w.trajectories) == 3)
            assert(len(storage_w.stores['snapshot0']) == 20)

            # this will not store since traj is already stored
            storage_w.trajectories.save(traj)
            assert(len(storage_w.snapshots) == 20)
            assert(len(storage_w.trajectories) == 3)
            assert(len(storage_w.stores['snapshot0']) == 20)

            # we saved in this order [0f, 8r, 6f, 7f, 9f, 5r, 4r, 3r, 2r, 1r ]
            # these are indices      [ 0, 17, 12, 14, 18,  3,  5,  7,  9, 11 ]

            storage_w.close()

            if use_cache:
                storage_r = paths.AnalysisStorage(fname)
            else:
                storage_r = paths.Storage(fname, 'r')
                storage_r.snapshots.set_caching(False)
                storage_r.stores['snapshot0'].set_caching(False)

            # check if the loaded trajectory is reproduced
            for s1, s2 in zip(traj, storage_r.trajectories[2]):
                compare_snapshot(s1, s2, True)

            # this is the expected order in which it is saved
            eff_traj = [
                traj[0],
                traj[8].reversed,
                traj[6],
                traj[7],
                traj[9],
                traj[5].reversed,
                traj[4].reversed,
                traj[3].reversed,
                traj[2].reversed,
                traj[1].reversed,
            ]

            # load from hidden and see, if the hidden store looks as expected
            # we open every second snapshot from the hidden store because the
            # ones in between correspond to the reversed ones

            hidden_snapshots = storage_r.stores['snapshot0'][:]
            for idx in range(10):
                s1 = eff_traj[idx]
                s1r = s1.reversed
                s2 = hidden_snapshots[2 * idx]
                s2r = hidden_snapshots[2 * idx + 1]
                compare_snapshot(s1, s2, True)
                compare_snapshot(s1r, s2r, True)

            storage_r.close()

    def test_load_save(self):
        for use_uuid in [True]:
            store = Storage(filename=self.filename, mode='w')
            assert(os.path.isfile(self.filename))

            store.save(self.template_snapshot)
            store.close()

            store = Storage(filename=self.filename, mode='a')
            loaded_template = store.snapshots[0]
            loaded_r = store.snapshots[1]

            compare_snapshot(loaded_template, self.template_snapshot, True)
            compare_snapshot(
                loaded_template.reversed,
                self.template_snapshot.reversed, True)
            compare_snapshot(loaded_r, self.template_snapshot.reversed)

            store.close()

    def test_proxy(self):
        for use_uuid in [True]:
            store = Storage(filename=self.filename, mode='w')
            assert(os.path.isfile(self.filename))

            tm = self.template_snapshot

            store.save(tm)

            px = store.snapshots.proxy(store.snapshots.index.list[0])

            # make sure that the proxy and
            assert(hash(px) == hash(tm))
            assert(px == tm)

            store.snapshots.cache.clear()
            s0 = store.snapshots[0]

            assert(hash(px) == hash(s0))
            assert(px == s0)

            compare_snapshot(px, tm)
            compare_snapshot(s0, tm)

            px = store.snapshots.proxy(store.snapshots.index.list[0])

            # make sure that after reloading it still works
            assert(hash(px) == hash(tm))
            assert(px == tm)

            store.close()

            store = Storage(filename=self.filename, mode='a')

            s1 = store.snapshots[0]

            store.close()

            # when loading only for uuid based storages you get the same id
            assert((hash(px) == hash(s1)) is use_uuid)
            assert((px == s1) is use_uuid)

    def test_mention_only(self):
        storage_w = paths.Storage(self.filename, "w")

        template = self.template_snapshot

        storage_w.snapshots.add_type(template)

        test_snap = self.traj[2]

        # only touch a new snapshot
        storage_w.snapshots.only_mention = True
        storage_w.snapshots.save(test_snap)

        # check that the snapshot is there
        assert(len(storage_w.snapshots) == 4)
        # in the memory uuid index
        assert(test_snap.__uuid__ in storage_w.snapshots.index)
        # and stored
        assert(test_snap.__uuid__ == storage_w.snapshots.vars['uuid'][1])

        # but no real snapshot has been stored
        # print len(storage_w.objects['snapshot0'])
        assert(len(storage_w.objects['snapshot0']) == 2)

        # switch on normal saving
        storage_w.snapshots.only_mention = False

        test_snap = self.traj[4]
        storage_w.snapshots.mention(test_snap)

        # check that the snapshot is there
        assert(len(storage_w.snapshots) == 6)
        # in the memory uuid index
        assert(test_snap.__uuid__ in storage_w.snapshots.index)
        # and stored
        assert(test_snap.__uuid__ == storage_w.snapshots.vars['uuid'][2])

        # but no real snapshot has been stored
        assert(len(storage_w.objects['snapshot0']) == 2)

        # try to now add it
        storage_w.snapshots.save(test_snap)

        # check that the snapshot is not stored again (only 3 snapshots)
        assert(len(storage_w.snapshots) == 6)
        assert(len(storage_w.objects['snapshot0']) == 4)

        # print storage_w.objects['snapshot0'][1].coordinates
        # print template.coordinates
        # print storage_w.objects['snapshot0'][0].coordinates
        # print test_snap.coordinates
        # print storage_w.objects['snapshot0'].vars['statics'][0].coordinates
        # print storage_w.objects['snapshot0'].vars['statics'][1].coordinates
        # print storage_w.objects['snapshot0'].index

        compare_snapshot(storage_w.objects['snapshot0'][4], test_snap)

        storage_w.close()

    def test_load_save_uuid(self):
        store = Storage(filename=self.filename, mode='w')
        assert(os.path.isfile(self.filename))

        store.save(self.template_snapshot)
        store.close()

        store = Storage(filename=self.filename, mode='a')
        loaded_template = store.snapshots[self.template_snapshot.__uuid__]
        loaded_r = store.snapshots[self.template_snapshot.reversed.__uuid__]

        compare_snapshot(loaded_template, self.template_snapshot, True)
        compare_snapshot(
            loaded_template.reversed,
            self.template_snapshot.reversed, True)
        compare_snapshot(loaded_r, self.template_snapshot.reversed)

        store.close()

    def test_load_save_toy(self):
        store = Storage(filename=self.filename, mode='w')
        assert(os.path.isfile(self.filename))

        store.save(self.toy_template)

        store.close()

        store = Storage(filename=self.filename, mode='a')
        loaded_template = store.snapshots[0]
        loaded_r = store.snapshots[1]

        compare_snapshot(loaded_template, self.toy_template, True)
        compare_snapshot(
            loaded_template.reversed,
            self.toy_template.reversed, True)
        compare_snapshot(loaded_r, self.toy_template.reversed)

        store.close()

    def test_reverse_bug(self):
        store = Storage(filename=self.filename,
                        mode='w')
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
        store = Storage(
            filename=self.filename, mode='w')

        assert(os.path.isfile(self.filename))
        assert(store.storage_version == paths.version.version)
