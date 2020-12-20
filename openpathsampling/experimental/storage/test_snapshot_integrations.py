# These are tests of the integration of various engines with the new storage
# system.

import pytest

import os

import numpy as np
import openpathsampling as paths

from openpathsampling.tests.test_helpers import data_filename

from openpathsampling.experimental.storage.ops_storage import (
    Storage, ops_class_info, ops_schema
)
from openpathsampling.experimental.simstore import SQLStorageBackend

class TestSnapshotIntegration(object):
    # TODO: no use of self anywhere in here; this might become bare test
    # functions
    def _make_storage(self, mode):
        Storage._known_storages = {}
        backend = SQLStorageBackend("test.sql", mode=mode)
        storage = Storage.from_backend(
            backend=backend,
            schema=ops_schema,
            class_info=ops_class_info
        )
        return storage

    def _make_gromacs_snap(self):
        # gromacs-specific
        test_dir = data_filename("gromacs_engine")
        engine = paths.engines.gromacs.Engine(
            gro="conf.gro",
            mdp="md.mdp",
            top="topol.top",
            options={},
            base_dir=test_dir,
            prefix="proj"
        )
        snap_file = os.path.join(test_dir, "project_trr", "0000000.trr")

        snap = paths.engines.gromacs.ExternalMDSnapshot(
            file_name=snap_file,
            file_position=2,
            engine=engine
        )
        return snap

    def _make_openmm_snap(self):
        mm = pytest.importorskip('simtk.openmm')
        traj_file = data_filename('ala_small_traj.pdb')
        traj = paths.engines.openmm.tools.ops_load_trajectory(traj_file)
        snap = traj[0]
        return snap

    @pytest.mark.parametrize('integration', ['gromacs', 'openmm'])
    def test_integration(self, integration):
        make_snap = {
            'gromacs': self._make_gromacs_snap,
            'openmm': self._make_openmm_snap,
        }[integration]
        snap = make_snap()

        storage = self._make_storage(mode='w')
        assert not storage.backend.has_table('snapshot0')
        storage.save(snap)
        assert storage.backend.has_table('snapshot0')
        storage.close()

        storage = self._make_storage(mode='r')
        assert storage.backend.has_table('snapshot0')
        snap2 = storage.snapshots.tables['snapshot0'][0]

        np.testing.assert_array_equal(snap.xyz, snap2.xyz)


