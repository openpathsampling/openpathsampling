import pytest

from openpathsampling.exports.trajectories.trrtrajectorywriter import *
from openpathsampling.integration_tools import HAS_MDTRAJ
from openpathsampling.utils.storage_interfaces import (
    LocalFileStorageInterface, MemoryStorageInterface
)
import numpy.testing as npt

def test_trr_trajectory_writer(ad_trajectory, tmp_path):
    if not HAS_MDTRAJ:
        pytest.skip("mdtraj is not available")

    import mdtraj as md
    outfile = tmp_path / "test.trr"
    writer = TRRTrajectoryWriter()
    assert not outfile.exists()
    writer(ad_trajectory, outfile)
    assert outfile.exists()

    trr = md.formats.TRRTrajectoryFile(str(outfile), mode='r')
    xyz, time, step, box, lambd, vel, forces = trr._read(
        int(len(ad_trajectory)),
        atom_indices=None,
        get_velocities=True,
    )

    assert forces is None
    npt.assert_allclose(xyz, ad_trajectory.xyz)
    npt.assert_allclose(vel, ad_trajectory.velocities)
    npt.assert_allclose(box, ad_trajectory.box_vectors)
    npt.assert_allclose(lambd, np.zeros(len(ad_trajectory)))
    npt.assert_allclose(time, np.zeros(len(ad_trajectory)))


@pytest.mark.parametrize("storage_type", ["memory", "local"])
def test_trr_trajectory_writer_storage(ad_trajectory, storage_type,
                                       tmp_path):
    if not HAS_MDTRAJ:
        pytest.skip("mdtraj is not available")

    import mdtraj as md
    storage = {
        "memory": MemoryStorageInterface(),
        "local": LocalFileStorageInterface(tmp_path / "storage"),
    }[storage_type]
    writer = TRRTrajectoryWriter()

    writer.write(ad_trajectory, storage, "test.trr")
    assert "test.trr" in storage

    outfile = tmp_path / "test.trr"
    storage.load("test.trr", outfile)
    trr = md.formats.TRRTrajectoryFile(str(outfile), mode='r')
    xyz, time, step, box, lambd, vel, forces = trr._read(
        int(len(ad_trajectory)),
        atom_indices=None,
        get_velocities=True,
    )

    assert forces is None
    npt.assert_allclose(xyz, ad_trajectory.xyz)
    npt.assert_allclose(vel, ad_trajectory.velocities)
    npt.assert_allclose(box, ad_trajectory.box_vectors)
    npt.assert_allclose(lambd, np.zeros(len(ad_trajectory)))
    npt.assert_allclose(time, np.zeros(len(ad_trajectory)))


def test_trr_writer_ext():
    if not HAS_MDTRAJ:
        pytest.skip("mdtraj is not available")

    writer = TRRTrajectoryWriter()
    assert writer.ext == "trr"
