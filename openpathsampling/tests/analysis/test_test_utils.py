import pytest

import collections

import openpathsampling as paths

from openpathsampling.tests.test_helpers import make_1d_traj

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
    if ensemble is not None:
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

def _setup_one_way_forward(scheme, accepted):
    partial_traj = make_trajectory(0.2, lower_bound=-0.1).reversed
    # * ensemble 0 is always accepted because the trial trajectory is
    #   shorter than the input trajectory
    # * ensemble 2 is always rejected because the trial trajectory doesn't
    #   cross the interface (therefore doesn't satisfy the ensemble)
    ensemble_idx = {True: 0, False: 2}[accepted]
    ensemble = scheme.network.sampling_ensembles[ensemble_idx]
    return ensemble, partial_traj

def _setup_one_way_backward(scheme, accepted):
    ensemble = scheme.network.sampling_ensembles[2]
    # * accepted trajectory always accepted because trial trajectory is
    #   shorter than the input trajectory
    # * rejected trajectory always rejected because it is B->B (does not
    #   satisfy the ensemble)
    partial_traj = {True: make_trajectory(0.1, lower_bound=-0.1).reversed,
                    False: make_trajectory(1.0, lower_bound=0.3)}[accepted]
    return ensemble, partial_traj

@pytest.mark.parametrize('direction', ['forward', 'backward'])
@pytest.mark.parametrize('accepted', [True, False])
def test_one_way_shooting_move(scheme, direction, accepted):
    t1 = make_trajectory(0.5)
    t2 = make_trajectory(1.0)
    init_conds = scheme.initial_conditions_from_trajectories([t1, t2])
    shooting_idx = 4

    ensemble, partial_traj = {
        'forward': _setup_one_way_forward(scheme, accepted),
        'backward': _setup_one_way_backward(scheme, accepted)
    }[direction]

    initial_traj = init_conds[ensemble].trajectory

    len_expected = {
        'forward': shooting_idx + 1 + len(partial_traj),
        'backward': len(partial_traj) + len(initial_traj) - shooting_idx
    }[direction]

    trial_shooting_idx = {'forward': shooting_idx,
                          'backward': len(partial_traj)}[direction]

    movetype = {'forward': MockForwardShooting,
                'backward': MockBackwardShooting}[direction]

    move = movetype(shooting_index=shooting_idx,
                    partial_traj=partial_traj,
                    scheme=scheme,
                    ensemble=ensemble)

    change = move(init_conds)

    # check that this looks like a move change from shooting
    assert len(change.trials) == 1
    trial_trajectory = change.trials[0].trajectory
    shooting_snapshot = change.subchanges[0].details.shooting_snapshot

    # check correctness of the results
    assert change.accepted is accepted
    assert len(trial_trajectory) == len_expected
    assert shooting_snapshot is initial_traj[shooting_idx]
    assert shooting_snapshot in trial_trajectory
    assert shooting_snapshot is trial_trajectory[trial_shooting_idx]


@pytest.mark.parametrize('accepted', [True, False, None])
def test_shooting_move_force_accept(scheme, accepted):
    # If the a given trial can be either accepted or rejected, then the
    # resulting change depends on the `accept` parameter of the mock move --
    # True accepts, False rejects, and None uses the default internal math.
    # Test this by trying each 25 times on a trial with a 50% acceptance
    # probability.
    # only test this with the forward shooting mover, since the logic is
    # shared with backward
    init_traj = make_trajectory(1.0)
    partial_traj = make_1d_traj([0.95] * (len(init_traj)-2) + [1.05])
    init_conds = scheme.initial_conditions_from_trajectories(init_traj)
    ensemble = scheme.network.sampling_ensembles[2]
    shooting_idx = len(init_traj) - 2

    move = MockForwardShooting(shooting_index=shooting_idx,
                               partial_traj=partial_traj,
                               accepted=accepted,
                               scheme=scheme,
                               ensemble=ensemble)

    n_attempts = 25
    changes = [move(init_conds) for _ in range(n_attempts)]
    results = collections.Counter(change.accepted for change in changes)
    if accepted is not None:
        assert results[accepted] == n_attempts
        assert not accepted not in results
    else:
        assert results[True] > 0
        assert results[False] > 0
        assert results[True] + results[False] == n_attempts

@pytest.mark.parametrize('accepted', [True, False])
def test_repex_move(scheme, accepted):
    t1 = make_trajectory(0.4)
    t2 = make_trajectory(1.0)
    init_conds = scheme.initial_conditions_from_trajectories([t1, t2])
    assert init_conds[0].trajectory is t1 and init_conds[1].trajectory is t1
    assert init_conds[2].trajectory is t2

    if accepted:
        e2, e1 = scheme.network.sampling_ensembles[:2]
    else:
        e2, e1 = scheme.network.sampling_ensembles[1:]

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
