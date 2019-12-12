"""
@author David W.H. Swenson
"""
from __future__ import absolute_import

from builtins import zip
from builtins import object
from .test_helpers import data_filename, assert_close_unit, md

import pytest
from nose.plugins.skip import SkipTest

import numpy as np

import openpathsampling.collectivevariable as op
import openpathsampling.engines.openmm as peng
from openpathsampling.netcdfplus import NetCDFPlus

try:
    from msmbuilder.featurizer import AtomPairsFeaturizer
except ImportError:
    has_msmbuilder = False
else:
    has_msmbuilder = True

import openpathsampling as paths
from openpathsampling.tests.test_helpers import make_1d_traj
import os


class TestFunctionCV(object):
    def setup(self):
        if not md:
            raise SkipTest("mdtraj not installed")
        self.mdtraj = md.load(data_filename("ala_small_traj.pdb"))
        _ = pytest.importorskip('simtk.unit')
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

    def test_pickle_external_cv(self):
        template = make_1d_traj([0.0])[0]
        cv = paths.FunctionCV("x", lambda snap: snap.coordinates[0][0])
        storage = paths.Storage("myfile.nc", "w", template)
        storage.save(cv)
        storage.close()

    def test_pickle_cv_with_imports(self):
        template = make_1d_traj([0.0])[0]

        def test_cv_func(snap):
            import math
            return math.ceil(snap.coordinates[0][0])

        cv = paths.FunctionCV("y", test_cv_func)
        storage = paths.Storage("myfile.nc", "w", template)
        storage.save(cv)
        storage.close()

    def test_dihedral_op(self):
        """ Create a dihedral order parameter """
        psi_atoms = [6, 8, 14, 16]
        dihedral_op = op.MDTrajFunctionCV(
            "psi",
            md.compute_dihedrals,
            topology=self.topology,
            indices=[psi_atoms])

        md_dihed = md.compute_dihedrals(self.mdtraj, indices=[psi_atoms])
        my_dihed = dihedral_op(self.traj_topology)

        np.testing.assert_allclose(
            md_dihed.reshape(md_dihed.shape[:-1]),
            my_dihed, rtol=10 ** -6, atol=10 ** -10)

    def test_atom_pair_featurizer(self):
        """ Create an atom pair collectivevariable using MSMSBuilder3 """

        if not has_msmbuilder:
            raise SkipTest("MSMBuilder not installed")
        atom_pairs = [[0, 1], [10, 14]]
        atom_pair_op = op.MSMBFeaturizerCV(
            "atom_pairs",
            AtomPairsFeaturizer,
            topology=self.topology,
            pair_indices=atom_pairs)

        md_distances = md.compute_distances(self.mdtraj, atom_pairs)

        my_distances = atom_pair_op(self.traj_topology)

        np.testing.assert_allclose(md_distances, my_distances, rtol=10 ** -6,
                                   atol=10 ** -10)

    def test_return_parameters_from_template(self):

        if not has_msmbuilder:
            raise SkipTest("MSMBuilder not installed")
        atom_pairs = [[0, 1], [10, 14]]
        atom_pair_op = op.MSMBFeaturizerCV(
            "atom_pairs",
            AtomPairsFeaturizer,
            topology=self.topology,
            pair_indices=atom_pairs)

        # little trick. We just pretend the atom_pairs_op is a function we want
        # to use it cannot be stored though, but for from_template it is enough

        params = NetCDFPlus.get_value_parameters(
            atom_pair_op(self.traj_topology[0]))

        assert params['var_type'] == 'numpy.float32'
        assert params['simtk_unit'] is None
        assert params['dimensions'] == tuple([2])

    def test_storage_cv_function(self):
        import os

        # test all combinations of (1) with and without UUIDs,
        # (2) using partial yes, no all of these must work
        for allow_incomplete in (True, False):

            # print '=========================================================='
            # print 'PARTIAL', allow_incomplete
            # print '=========================================================='

            fname = data_filename("cv_storage_test.nc")
            if os.path.isfile(fname):
                os.remove(fname)

            traj = paths.Trajectory(list(self.traj_simple))
            template = traj[0]

            storage_w = paths.Storage(fname, "w")
            storage_w.snapshots.save(template)

            cv1 = paths.CoordinateFunctionCV(
                'f1',
                lambda x: x.coordinates[0]
            ).with_diskcache(
                allow_incomplete=allow_incomplete
            )

            storage_w.save(cv1)

            # let's mess up the order in which we save and
            # include reversed ones as well
            assert (len(storage_w.snapshots) == 2)
            storage_w.trajectories.save(traj[3:])
            assert (len(storage_w.snapshots) == 16)
            storage_w.snapshots.save(traj[1].reversed)
            assert (len(storage_w.snapshots) == 18)
            storage_w.trajectories.save(traj.reversed)
            assert (len(storage_w.snapshots) == 20)

            # this should be ignored for all is saved already
            storage_w.trajectories.save(traj)
            storage_w.close()

            storage_r = paths.AnalysisStorage(fname)
            rcv1 = storage_r.cvs['f1']

            assert(rcv1._store_dict)

            cv_cache = rcv1._store_dict.value_store

            assert (cv_cache.allow_incomplete == allow_incomplete)

            for idx, snap in enumerate(storage_r.trajectories[1]):
                # print idx, snap
                # if hasattr(snap, '_idx'):
                #     print 'Proxy IDX', snap._idx

                # print 'ITEMS', storage_r.snapshots.index.items()
                # print snap, type(snap), snap.__dict__

                # print snap.__uuid__
                # print snap.reversed.__uuid__
                # print snap.create_reversed().__uuid__
                #
                # print 'POS', cv_cache.object_pos(snap),
                # print 'POS', storage_r.snapshots.pos(snap),
                # print 'POS', storage_r.snapshots.index[snap]
                #
                # print 'POS', cv_cache.object_pos(snap.reversed),
                # print 'POS', storage_r.snapshots.pos(snap.reversed),
                # print 'POS', storage_r.snapshots.index[snap.reversed]

                # if len(cv_cache.cache._chunkdict) > 0:
                #
                #     if allow_incomplete:
                #         print cv_cache.index
                #         print cv_cache.vars['value'][:]
                #
                #     for n, v in enumerate(cv_cache.cache._chunkdict[0]):
                #         print n, v
                #
                # print cv1(snap)
                # print cv1(snap.reversed)
                # print cv_cache[snap]
                #
                # print cv_cache[snap.reversed]

                if not allow_incomplete or cv_cache[snap] is not None:
                    assert_close_unit(cv_cache[snap], cv1(snap))
                    assert_close_unit(
                        cv_cache[snap.reversed],
                        cv1(snap.reversed))

            storage_r.close()

            if os.path.isfile(fname):
                os.remove(fname)

    def test_storage_sync_and_complete(self):
        import os

        # test all combinations of (1) with and without UUIDs,
        # (2) using partial yes, no all of these must work

        allow_incomplete = True

        # print
        # print

        for use_uuid in [True]:

            # print '=========================================================='
            # print 'UUID', use_uuid
            # print '=========================================================='

            fname = data_filename("cv_storage_test.nc")
            if os.path.isfile(fname):
                os.remove(fname)

            traj = paths.Trajectory(list(self.traj_simple))
            template = traj[0]

            storage_w = paths.Storage(fname, "w")
            storage_w.snapshots.save(template)

            cv1 = paths.CoordinateFunctionCV(
                'f1',
                lambda snapshot: snapshot.coordinates[0]
            ).with_diskcache(
                allow_incomplete=allow_incomplete
            )

            # let's mess up the order in which we save and include
            # reversed ones as well
            assert (len(storage_w.snapshots) == 2)
            storage_w.trajectories.save(traj[3:])
            assert (len(storage_w.snapshots) == 16)
            storage_w.snapshots.save(traj[1].reversed)
            assert (len(storage_w.snapshots) == 18)
            storage_w.trajectories.save(traj.reversed)
            assert (len(storage_w.snapshots) == 20)

            storage_w.save(cv1)

            store = storage_w.cvs.cache_store(cv1)
            assert (len(store.vars['value']) == 0)

            storage_w.snapshots.complete_cv(cv1)
            assert (len(store.vars['value']) == 10)

            # check if stored values match computed ones
            for idx, value in zip(
                    store.variables['index'][:],
                    store.vars['value']):
                snap = storage_w.snapshots[
                    storage_w.snapshots.vars['uuid'][idx]]

                # print(snap, snap.__uuid__, value)
                assert_close_unit(cv1(snap), value)

            storage_w.close()

            if os.path.isfile(fname):
                os.remove(fname)

    def test_storage_sync(self):
        import os

        # test all combinations of (1) with and without UUIDs,
        # (2) using partial yes; all of these must work

        allow_incomplete = True

        # print
        # print

        for use_uuid in [True]:

            # print '=========================================================='
            # print 'UUID', use_uuid
            # print '=========================================================='

            fname = data_filename("cv_storage_test.nc")
            if os.path.isfile(fname):
                os.remove(fname)

            traj = paths.Trajectory(list(self.traj_simple))
            template = traj[0]

            storage_w = paths.Storage(fname, "w")
            storage_w.snapshots.save(template)

            cv1 = paths.CoordinateFunctionCV(
                'f1',
                lambda snapshot: snapshot.coordinates[0]
            ).with_diskcache(
                allow_incomplete=allow_incomplete)

            # let's mess up the order in which we save and
            # include reversed ones as well
            assert (len(storage_w.snapshots) == 2)
            storage_w.trajectories.save(traj[6:])
            assert (len(storage_w.snapshots) == 10)
            storage_w.snapshots.save(traj[1].reversed)
            assert (len(storage_w.snapshots) == 12)

            storage_w.save(cv1)

            store = storage_w.cvs.cache_store(cv1)
            assert (len(store.vars['value']) == 0)
            storage_w.snapshots.sync_cv(cv1)

            # nothing added to the cache so no changes
            assert (len(store.vars['value']) == 1)

            # fill the cache
            _ = cv1(traj)

            storage_w.snapshots.sync_cv(cv1)

            # should match the number of stored snapshots
            assert (len(store.vars['value']) == 6)

            # save the rest
            storage_w.trajectories.save(traj.reversed)
            assert (len(storage_w.snapshots) == 20)

            # should still be unchanged
            assert (len(store.vars['value']) == 6)

            # this should store the remaining CV values

            storage_w.snapshots.sync_cv(cv1)
            assert (len(store.vars['value']) == 10)

            # check if the values match
            for idx, value in zip(
                    store.variables['index'][:],
                    store.vars['value']):
                snap = storage_w.snapshots[
                    storage_w.snapshots.vars['uuid'][idx]]

                assert_close_unit(cv1(snap), value)

            storage_w.close()

            if os.path.isfile(fname):
                os.remove(fname)
