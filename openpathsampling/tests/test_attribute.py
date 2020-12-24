"""
@author David W.H. Swenson
"""

from .test_helpers import data_filename, assert_close_unit, make_1d_traj, md
import pytest

from nose.plugins.skip import SkipTest

import openpathsampling.engines.openmm as peng
from openpathsampling.netcdfplus import FunctionPseudoAttribute

import openpathsampling as paths
import os


class TestFunctionPseudoAttribute(object):
    def setup(self):
        if not md:
            raise SkipTest("mdtraj not installed")
        self.mdtraj = md.load(data_filename("ala_small_traj.pdb"))
        pytest.importorskip("simtk.unit")
        self.traj_topology = peng.trajectory_from_mdtraj(self.mdtraj)
        self.traj_simple = peng.trajectory_from_mdtraj(
            self.mdtraj,
            simple_topology=True)

        self.topology = self.traj_topology[0].engine.topology

        if os.path.isfile("myfile.nc"):
            os.remove("myfile.nc")

    def teardown(self):
        if os.path.isfile("myfile.nc"):
            os.remove("myfile.nc")

    def test_pickle_external_attr(self):
        template = make_1d_traj([0.0])[0]
        attr = FunctionPseudoAttribute("x", paths.Trajectory, lambda x: x)
        storage = paths.Storage("myfile.nc", "w", template)
        storage.save(attr)
        storage.close()

    def test_storage_attribute_function(self):
        import os

        # test all combinations of (1) with and without UUIDs,
        # (2) using partial yes, no all of these must work
        for allow_incomplete in (True, False):

            # print('ALLOW INCOMPLETE', allow_incomplete)

            fname = data_filename("attr_storage_test.nc")
            if os.path.isfile(fname):
                os.remove(fname)

            traj = paths.Trajectory(list(self.traj_simple))
            template = traj[0]

            storage_w = paths.Storage(fname, "w")
            storage_w.snapshots.save(template)
            storage_w.trajectories.save(traj)

            # compute distance in x[0]
            attr1 = FunctionPseudoAttribute(
                'f1',
                paths.Trajectory,
                lambda x: x[0].coordinates[0] - x[-1].coordinates[0]
            ).with_diskcache(
                allow_incomplete=allow_incomplete
            )

            storage_w.save(attr1)

            # let's mess up the order in which we save and
            # include reversed ones as well
            storage_w.trajectories.save(traj[3:])
            storage_w.snapshots.save(traj[1].reversed)
            storage_w.trajectories.save(traj.reversed)

            # this should be ignored for all is saved already
            storage_w.trajectories.save(traj)
            storage_w.close()

            storage_r = paths.Storage(fname, 'r')
            rattr1 = storage_r.attributes['f1']

            assert(rattr1._store_dict)

            attr_cache = rattr1._store_dict.value_store

            assert (attr_cache.allow_incomplete == allow_incomplete)

            for idx, traj in enumerate(storage_r.trajectories):
                if not allow_incomplete or attr_cache[traj] is not None:
                    assert_close_unit(attr_cache[traj], attr1(traj))

            storage_r.close()

            if os.path.isfile(fname):
                os.remove(fname)

    def test_storage_sync_and_complete(self):
        import os

        # test all combinations of (1) with and without UUIDs,
        # (2) using partial yes, no all of these must work

        allow_incomplete = True

        fname = data_filename("attr_storage_test.nc")
        if os.path.isfile(fname):
            os.remove(fname)

        traj = paths.Trajectory(list(self.traj_simple))
        template = traj[0]

        storage_w = paths.Storage(fname, "w")
        storage_w.snapshots.save(template)
        storage_w.trajectories.save(traj)

        # compute distance in x[0]
        attr1 = FunctionPseudoAttribute(
            'f1',
            paths.Trajectory,
            lambda x: x[0].coordinates[0] - x[-1].coordinates[0]
        ).with_diskcache(
            allow_incomplete=allow_incomplete
        )

        storage_w.save(attr1)

        store = storage_w.attributes.cache_store(attr1)
        storage_w.trajectories.complete_attribute(attr1)

        # check if stored values match computed ones
        for idx, value in zip(
                store.variables['index'][:],
                store.vars['value']):
            traj = storage_w.trajectories[
                storage_w.trajectories.vars['uuid'][idx]]

            assert_close_unit(attr1(traj), value)

        storage_w.close()

        if os.path.isfile(fname):
            os.remove(fname)

    def test_storage_sync(self):
        import os

        # test all combinations of (1) with and without UUIDs,
        # (2) using partial yes; all of these must work

        allow_incomplete = True

        fname = data_filename("attr_storage_test.nc")
        if os.path.isfile(fname):
            os.remove(fname)

        traj = paths.Trajectory(list(self.traj_simple))
        template = traj[0]

        storage_w = paths.Storage(fname, "w")
        storage_w.snapshots.save(template)
        storage_w.trajectories.save(traj)

        # compute distance in x[0]
        attr1 = FunctionPseudoAttribute(
            'f1',
            paths.Trajectory,
            lambda x: x[0].coordinates[0] - x[-1].coordinates[0]
        ).with_diskcache(
            allow_incomplete=allow_incomplete
        )

        storage_w.save(attr1)

        store = storage_w.attributes.cache_store(attr1)
        storage_w.trajectories.sync_attribute(attr1)

        # fill the cache
        _ = attr1(traj)

        storage_w.trajectories.sync_attribute(attr1)

        # should match the number of stored snapshots
        # assert (len(store.vars['value']) == 6)

        # save the rest
        storage_w.trajectories.save(traj.reversed)
        # assert (len(storage_w.snapshots) == 20)

        # should still be unchanged
        # assert (len(store.vars['value']) == 6)

        # this should store the remaining attr values

        storage_w.trajectories.sync_attribute(attr1)
        # assert (len(store.vars['value']) == 10)

        # check if the values match
        for idx, value in zip(
                store.variables['index'][:],
                store.vars['value']):
            snap = storage_w.trajectories[
                storage_w.trajectories.vars['uuid'][idx]]

            assert_close_unit(attr1(snap), value)

        storage_w.close()

        if os.path.isfile(fname):
            os.remove(fname)
