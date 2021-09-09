import pytest

import openpathsampling as paths

from openpathsampling.tests.analysis.utils import *

@pytest.mark.parametrize('maxval', [0.5, 1.0])
def test_make_trajectory(maxval):
    traj = make_trajectory(maxval)
    xvals = traj.xyz[:,0,0]
    assert maxval <= max(xvals) < maxval + 0.1
    assert -0.1 <= min(xvals) < 0.0
    assert -0.1 <= xvals[0] < 0.0
    if maxval == 0.5:
        assert len(traj) == 13
        assert -0.1 <= xvals[-1] < 0.0
    elif maxval == 1.0:
        assert len(traj) == 12
        assert 1.0 <= xvals[-1] < 1.1
    else:
        raise RuntimeError("This shouldn't happen")

@pytest.mark.parametrize('ensemble', [None, 'ensemble'])
def test_select_by_input_ensembles(scheme, ensemble):
    pass

def test_wrap_one_way_shooting_mover(scheme):
    pass

def test_wrap_repex_move():
    pass
