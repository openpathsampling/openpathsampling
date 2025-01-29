import pytest

from openpathsampling.exports.trajectories.mdtrajtrajectorywriter import *
from openpathsampling.integration_tools import HAS_MDTRAJ
import numpy.testing as npt

def test_mdtraj_trajectory_writer(ad_trajectory, ad_grofile, tmp_path):
    if not HAS_MDTRAJ:
        pytest.skip("mdtraj is not available")

    import mdtraj as md
    outfile = tmp_path / "test.xtc"
    writer = MDTrajTrajectoryWriter()
    assert not outfile.exists()
    writer(ad_trajectory, outfile)
    assert outfile.exists()

    reloaded = md.load(str(outfile), top=ad_grofile)
    assert len(reloaded) == len(ad_trajectory)
    assert reloaded.xyz.shape == ad_trajectory.xyz.shape
    npt.assert_allclose(reloaded.xyz, ad_trajectory.xyz, atol=1e-3)


@pytest.mark.parametrize("selection", [
    "backbone", 
    [4,  5,  6,  8, 14, 15, 16, 18]
])
def test_subtrajectory_selection(ad_trajectory, selection, tmp_path):
    if not HAS_MDTRAJ:
        pytest.skip("mdtraj is not available")

    import mdtraj as md
    outfile = tmp_path / "test.xtc"
    writer = MDTrajTrajectoryWriter(mdtraj_selection=selection)
    assert not outfile.exists()
    writer(ad_trajectory, outfile)
    assert outfile.exists()

    # save file for the modified topology
    mdt = ad_trajectory.to_mdtraj()
    if isinstance(selection, str):
        subset = mdt.topology.select(selection)
    else:
        subset = selection

    subtraj = mdt.atom_slice(subset)
    subgro = tmp_path / "subset.gro"
    subtraj[0].save(str(subgro))

    reloaded = md.load(str(outfile), top=subgro)
    assert len(reloaded) == len(ad_trajectory)
    assert reloaded.xyz.shape == subtraj.xyz.shape
    assert reloaded.xyz.shape != ad_trajectory.xyz.shape
    npt.assert_allclose(reloaded.xyz, subtraj.xyz, atol=1e-3)


def test_mdtraj_trajectory_writer_selection_error(ad_trajectory, tmp_path):
    if not HAS_MDTRAJ:
        pytest.skip("mdtraj is not available")

    writer = MDTrajTrajectoryWriter(mdtraj_selection=object())
    with pytest.raises(TypeError):
        writer(ad_trajectory, tmp_path / "test.xtc")
