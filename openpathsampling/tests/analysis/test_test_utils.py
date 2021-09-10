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

@pytest.mark.parametrize('maxval', [0.5, 1.0])
def test_make_trajectory_lower_bound(maxval):
    traj = make_trajectory(maxval, 0.3)
    xvals = traj.xyz[:,0,0]
    assert 0.3 < min(xvals) <= 0.4
    assert max(xvals) == xvals[-1]
    if maxval == 0.5:
        assert len(traj) == 3
        assert 0.5 < max(xvals) <= 0.6
    elif maxval == 1.0:
        assert len(traj) == 8
        assert 1.0 < max(xvals) <= 1.1
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

def test_random_choice_mover(scheme):
    mover = scheme.movers['shooting'][0]
    target_change = Mock(mover=mover)
    shooting_choosers = [group for group in scheme.root_mover.submovers
                         if mover in group.submovers]
    assert len(shooting_choosers) == 1
    shooting_chooser = shooting_choosers[0]
    # make sure we got the right one
    assert shooting_chooser.name == "ShootingChooser"

    random_choice = MockRandomChoiceMover(shooting_chooser, target_change)
    change = random_choice('foo')  # exact inputs are ignored
    assert change.mover is shooting_chooser
    assert change.subchanges[0] is target_change

@pytest.mark.parametrize('accepted', [True, False])
def test_forward_shooting_move(scheme, accepted):
    pytest.skip()

@pytest.mark.parametrize('accepted', [True, False])
def test_backward_shooting_move(scheme, accepted):
    pytest.skip()

@pytest.mark.parametrize('accepted', [True, False])
def test_repex_move(scheme, accepted):
    t1 = make_trajectory(0.4)
    t2 = make_trajectory(1.0)
    init_conds = scheme.initial_conditions_from_trajectories([t1, t2])
    assert init_conds[0].trajectory is t1 and init_conds[1].trajectory is t1
    assert init_conds[2].trajectory is t2

    if accepted:
        e1, e2 = scheme.network.sampling_ensembles[:2]
    else:
        e1, e2 = scheme.network.sampling_ensembles[1:]

    ensembles = e1, e2
    move = MockRepex(scheme, ensembles)
    change = move(init_conds)
    assert change.accepted is accepted
    assert isinstance(change.canonical.mover, paths.ReplicaExchangeMover)

@pytest.mark.parametrize('accepted', [True, False])
def test_mock_pathreversal(scheme, accepted):
    move = MockPathReversal(scheme)
    maxval = {True: 0.6, False: 1.0}[accepted]
    init_traj = make_trajectory(maxval)
    init_conds = scheme.initial_conditions_from_trajectories(init_traj)
    change = move(init_conds)
    assert change.accepted is accepted
    assert isinstance(change.canonical.mover, paths.PathReversalMover)

@pytest.mark.parametrize('accepted', [True, False])
def test_wrap_org_by_group(scheme, accepted):
    # the specific mover doesn't matter; path reversal is easy
    move = MockPathReversal(scheme)
    maxval = {True: 0.6, False: 1.0}[accepted]
    init_traj = make_trajectory(maxval)
    init_conds = scheme.initial_conditions_from_trajectories(init_traj)
    canonical = move(init_conds)
    wrapped = move.wrap_org_by_group(canonical, init_conds)
    assert wrapped.mover is scheme.root_mover
    assert wrapped.subchanges[0].mover.name == "PathreversalChooser"
    assert wrapped.canonical is canonical
    assert wrapped.subchanges[0].subchanges[0] is canonical

def test_steps():
    pytest.skip()