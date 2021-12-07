import pytest

from openpathsampling.tests.analysis.utils.fixture_classes import *


@pytest.mark.parametrize('bounds', [(-1, 10), (-1, 5)])
def test_make_trajectory(bounds, default_two_state_tps):
    make_trajectory = default_two_state_tps.make_trajectory
    lower, upper = bounds
    traj = make_trajectory(lower, upper)
    xvals = traj.xyz[:, 0, 0]
    assert upper <= max(xvals) < upper + 1
    assert lower <= min(xvals) < lower + 1
    expected_len = upper - lower + 1
    assert len(traj) == expected_len


@pytest.mark.parametrize('lower_bound', [None, -2])
def test_make_tis_trajectory(lower_bound, default_unidirectional_tis):
    make_tis_trajectory = default_unidirectional_tis.make_tis_trajectory
    cv_max = 5
    kwargs = {'lower_bound': lower_bound} if lower_bound is not None else {}
    traj = make_tis_trajectory(cv_max, **kwargs)
    xvals = traj.xyz[:, 0, 0]

    if lower_bound is None:
        lower_bound = -1

    assert lower_bound <= min(xvals) < lower_bound + 1
    assert cv_max <= max(xvals) < cv_max + 1
    expected_len = 2 * (cv_max - lower_bound) + 1
    assert len(traj) == expected_len


def test_make_tis_trajectory_transition(default_unidirectional_tis):
    make_tis_trajectory = default_unidirectional_tis.make_tis_trajectory
    traj = make_tis_trajectory(10)
    xvals = traj.xyz[:, 0, 0]
    assert max(xvals) == xvals[-1]
