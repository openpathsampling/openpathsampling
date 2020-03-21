import pytest
from numpy import testing as npt
import openpathsampling as paths
import os

try:
    import openmmtools
except ImportError:
    openmmtools = None


class TestStaticContainerStore(object):
    def teardown(self):
        if os.path.isfile("test.nc"):
            os.remove("test.nc")

    def test_store_nonperiodic(self):
        if not openmmtools:
            pytest.skip("Requires OpenMMTools for testing")
        testsystem = openmmtools.testsystems.AlanineDipeptideVacuum()
        snap = paths.engines.openmm.snapshot_from_testsystem(testsystem,
                                                             periodic=False)
        storage = paths.Storage("test.nc", 'w')
        storage.save(snap)
        storage.close()
        load = paths.Storage("test.nc", 'r')
        reloaded = load.snapshots[0]
        npt.assert_array_equal(snap.coordinates, reloaded.coordinates)
        npt.assert_array_equal(snap.box_vectors, reloaded.box_vectors)
        assert snap.box_vectors is None
        assert reloaded.box_vectors is None
