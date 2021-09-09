import pytest

import openpathsampling as paths

from openpathsampling.tests.analysis.utils import *
from openpathsampling.tests.analysis.utils import _select_by_input_ensembles

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
    movers = scheme.movers['shooting']
    if ensemble == 'ensemble':
        ensemble = movers[0].input_ensembles[0]

    selected = _select_by_input_ensembles(movers, ensemble)
    assert selected in movers  # this is enough for ensemble=None
    if ensemble == 'ensemble':
        assert selected == movers[0]

def test_forward_shooting_move(scheme):
    pytest.skip()

def test_backward_shooting_move(scheme):
    pytest.skip()

@pytest.mark.parametrize('accepted', [True, False])
def test_repex_move(scheme, accepted):
    pytest.skip()

@pytest.mark.parametrize('accepted', [True, False])
def test_mock_pathreversal(scheme, accepted):
    move = MockPathReversal(scheme)
    maxval = {True: 0.6, False: 1.0}[accepted]
    init_traj = make_trajectory(maxval)
    init_conds = scheme.initial_conditions_from_trajectories(init_traj)
    change = move(init_conds)
    assert change.accepted is accepted
    assert isinstance(change.canonical.mover, paths.PathReversalMover)

def test_steps():
    pytest.skip()
