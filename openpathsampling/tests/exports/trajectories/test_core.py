from openpathsampling.exports.trajectories import *

from openpathsampling.tests.test_helpers import (
    data_filename, make_1d_traj
)
import openpathsampling as paths

import pytest
import pathlib


@pytest.fixture
def _ad_gmx_engine():
    engine = paths.engines.gromacs.Engine(
        gro="conf.gro",
        mdp="md.mdp",
        top="topol.top",
        options = {
            'mdrun_args': '-nt 1',
            'grompp_args': '-maxwarn 2',
        },
        base_dir=data_filename("gromacs_engine"),
        prefix="project"
    )
    return engine


@pytest.fixture
def ad_trajpath():
    test_dir = pathlib.Path(data_filename("gromacs_engine"))
    trajfile = test_dir / "project_trr/0000000.trr"
    return trajfile


@pytest.fixture
def ad_trajectory(ad_trajpath, _ad_gmx_engine):
    traj = paths.Trajectory([
        _ad_gmx_engine.read_frame_from_file(str(ad_trajpath), i)
        for i in [0, 1, 2]
    ])
    return traj


@pytest.fixture
def toy_trajectory():
    return make_1d_traj([1.0, 2.0, 3.0])


class TrajectoryWriterTestBase:
    def _read_trajectory(self, filename):
        raise NotImplementedError()

    def _test_trajectory(self, trajectory, outfile):
        assert not outfile.exists()
        self.writer(trajectory, outfile)
        assert outfile.exists()

        traj = self._read_trajectory(outfile)
        assert len(traj) == len(trajectory)
        assert traj == trajectory

    def _test_call_outfile_exists(self, trajectory, outfile):
        outfile.touch()
        with pytest.raises(FileExistsError):
            self.writer(trajectory, outfile)

    def _test_call_outfile_exists_force(self, trajectory, outfile):
        outfile.touch()
        assert outfile.exists()
        assert len(outfile.read_bytes()) == 0
        self.writer(trajectory, outfile, force=True)
        assert outfile.exists()
        assert len(outfile.read_bytes()) > 0

    def test_call(self, trajectory_fixture, request, tmp_path):
        raise NotImplementedError()

    def test_call_outfile_exists(self, request, tmp_path):
        raise NotImplementedError()

    def test_call_outfile_exists_force(self, request, tmp_path):
        raise NotImplementedError()


class TestSimStoreTrajectoryWriter(TrajectoryWriterTestBase):
    def setup_method(self):
        self.writer = SimStoreTrajectoryWriter()

    def _read_trajectory(self, filename):
        from openpathsampling.experimental.storage import Storage
        storage = Storage(filename, mode='r')
        return storage.trajectories[0]

    @pytest.mark.parametrize("trajectory_fixture", [
        "ad_trajectory",
        "toy_trajectory",
    ])
    def test_call(self, trajectory_fixture, request, tmp_path):
        # monkey patch for SimStore
        trajectory = request.getfixturevalue(trajectory_fixture)
        import openpathsampling as paths
        from openpathsampling.experimental.storage import monkey_patches
        paths = monkey_patches.monkey_patch_all(paths)

        outfile = tmp_path / f"{trajectory_fixture}.nc"
        self._test_trajectory(trajectory, outfile)

        # undo the monkey patch
        monkey_patches.unpatch(paths)

    def test_call_outfile_exists(self, request, tmp_path):
        trajectory = request.getfixturevalue("ad_trajectory")
        self._test_call_outfile_exists(trajectory, tmp_path / "test.db")

    def test_call_outfile_exists_force(self, request, tmp_path):
        trajectory = request.getfixturevalue("ad_trajectory")
        self._test_call_outfile_exists_force(trajectory, tmp_path / "test.db")
