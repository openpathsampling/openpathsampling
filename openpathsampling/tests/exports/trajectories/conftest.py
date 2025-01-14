import pytest
import pathlib
from openpathsampling.tests.test_helpers import data_filename
import openpathsampling as paths

@pytest.fixture
def ad_trajpath():
    test_dir = pathlib.Path(data_filename("gromacs_engine"))
    trajfile = test_dir / "project_trr/0000000.trr"
    return trajfile

@pytest.fixture
def ad_grofile():
    test_dir = pathlib.Path(data_filename("gromacs_engine"))
    topfile = test_dir / "conf.gro"
    return str(topfile)


@pytest.fixture
def ad_trajectory(ad_trajpath):
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

    traj = paths.Trajectory([
        engine.read_frame_from_file(str(ad_trajpath), i)
        for i in [0, 1, 2]
    ])
    return traj
