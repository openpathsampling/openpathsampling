import pytest

try:
    from unittest.mock import Mock
except ImportError:  # pragma: no cover
    from mock import Mock

import collections
import openpathsampling as paths

from openpathsampling.tests.test_helpers import make_1d_traj

from openpathsampling.tests.analysis.utils.fixture_classes import (
    TPSSystemFixture
)
from openpathsampling.tests.analysis.utils.mock_movers import *
from openpathsampling.tests.analysis.utils.mock_movers import (
    _select_by_input_ensembles, _get_only
)


def test_get_only():
    test_list = ['foo', 'bar', 'baz', 'qux']
    only = _get_only(iterable=test_list,
                     condition=lambda s: s.startswith('f'),
                     error_msg="string starting with f")
    assert only == 'foo'


def test_get_only_error():
    test_list = ['foo', 'bar', 'baz', 'qux']
    with pytest.raises(AnalysisTestSetupError, match="found 2"):
        _get_only(test_list,
                  condition=lambda s: s.startswith('b'),
                  error_msg="string starting with b")

    with pytest.raises(AnalysisTestSetupError, match="found 0"):
        _get_only(test_list,
                  condition=lambda s: s.startswith('z'),
                  error_msg="string starting with z")


@pytest.mark.parametrize('ensemble', [None, 'ensemble'])
def test_select_by_input_ensembles(default_unidirectional_tis, ensemble):
    scheme = default_unidirectional_tis.scheme
    movers = scheme.movers['shooting']
    if ensemble == 'ensemble':
        ensemble = movers[0].input_ensembles[0]

    selected = _select_by_input_ensembles(movers, ensemble)
    assert selected in movers  # this is enough for ensemble=None
    if ensemble is not None:
        assert selected == movers[0]


def test_random_choice_mover(default_unidirectional_tis):
    scheme = default_unidirectional_tis.scheme
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


def test_random_choice_mover_error(default_unidirectional_tis):
    scheme = default_unidirectional_tis.scheme
    mover = scheme.movers['repex'][0]
    target_change = Mock(mover=mover)
    shoot_mover = scheme.movers['shooting'][0]
    shooting_choosers = [group for group in scheme.root_mover.submovers
                         if shoot_mover in group.submovers]
    assert len(shooting_choosers) == 1
    shooting_chooser = shooting_choosers[0]
    # make sure we got the right one
    assert shooting_chooser.name == "ShootingChooser"

    random_choice = MockRandomChoiceMover(shooting_chooser, target_change)
    with pytest.raises(AnalysisTestSetupError, match="not found as sub"):
        random_choice('foo')


def _setup_one_way_forward(scheme, accepted):
    make_trajectory = TPSSystemFixture.make_trajectory
    partial_traj = make_trajectory(-1, 2).reversed
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
    make_trajectory = TPSSystemFixture.make_trajectory
    partial_traj = {True: make_trajectory(-1, 1).reversed,
                    False: make_trajectory(-3, 10)}[accepted]
    return ensemble, partial_traj


@pytest.mark.parametrize('direction', ['forward', 'backward'])
@pytest.mark.parametrize('accepted', [True, False])
def test_one_way_shooting_move(default_unidirectional_tis, direction,
                               accepted):
    scheme = default_unidirectional_tis.scheme
    t1 = default_unidirectional_tis.make_tis_trajectory(5)
    t2 = default_unidirectional_tis.make_tis_trajectory(10)
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

    assert move.ensemble == ensemble

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
def test_shooting_move_force_accept(default_unidirectional_tis, accepted):
    # If the a given trial can be either accepted or rejected, then the
    # resulting change depends on the `accept` parameter of the mock move --
    # True accepts, False rejects, and None uses the default internal math.
    # Test this by trying each 25 times on a trial with a 50% acceptance
    # probability.
    # only test this with the forward shooting mover, since the logic is
    # shared with backward
    scheme = default_unidirectional_tis.scheme
    init_traj = default_unidirectional_tis.make_tis_trajectory(10)
    partial_traj = make_1d_traj([9.5] * (len(init_traj)-2) + [10.5])
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
    # TODO: add assert to check that the pick probability is as expected?
    results = collections.Counter(change.accepted for change in changes)
    if accepted is not None:
        assert results[accepted] == n_attempts
        assert not accepted not in results
    else:
        assert results[True] > 0
        assert results[False] > 0
        assert results[True] + results[False] == n_attempts


def test_reject_nonsense_forced_acceptance(default_unidirectional_tis):
    # if a user requires that the shooting move be accepted, but the given
    # trajectory cannot be accepted because it doesn't match the ensemble,
    # we should raise an error
    scheme = default_unidirectional_tis.scheme
    init_traj = default_unidirectional_tis.make_tis_trajectory(10)
    init_conds = scheme.initial_conditions_from_trajectories(init_traj)
    # get a setup that should never be accepted
    ensemble, partial_traj = _setup_one_way_backward(scheme, accepted=False)
    shooting_idx = 4
    move = MockBackwardShooting(shooting_index=shooting_idx,
                                partial_traj=partial_traj,
                                accepted=True,
                                scheme=scheme,
                                ensemble=ensemble)

    with pytest.raises(AnalysisTestSetupError, match="force acceptance"):
        move(init_conds)


@pytest.mark.parametrize('accepted', [True, False])
def test_repex_move(default_unidirectional_tis, accepted):
    scheme = default_unidirectional_tis.scheme
    t1 = default_unidirectional_tis.make_tis_trajectory(4)
    t2 = default_unidirectional_tis.make_tis_trajectory(10)
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
def test_mock_pathreversal(default_unidirectional_tis, accepted):
    scheme = default_unidirectional_tis.scheme
    move = MockPathReversal(scheme)
    maxval = {True: 6, False: 10}[accepted]
    init_traj = default_unidirectional_tis.make_tis_trajectory(maxval)
    init_conds = scheme.initial_conditions_from_trajectories(init_traj)
    change = move(init_conds)
    assert change.accepted is accepted
    assert isinstance(change.canonical.mover, paths.PathReversalMover)


@pytest.mark.parametrize('accepted', [True, False])
def test_wrap_org_by_group(default_unidirectional_tis, accepted):
    # the specific mover doesn't matter; path reversal is easy
    scheme = default_unidirectional_tis.scheme
    move = MockPathReversal(scheme)
    maxval = {True: 6, False: 10}[accepted]
    init_traj = default_unidirectional_tis.make_tis_trajectory(maxval)
    init_conds = scheme.initial_conditions_from_trajectories(init_traj)
    canonical = move(init_conds)
    wrapped = move.wrap_org_by_group(canonical, init_conds)
    assert wrapped.mover is scheme.root_mover
    assert wrapped.subchanges[0].mover.name == "PathreversalChooser"
    assert wrapped.canonical is canonical
    assert wrapped.subchanges[0].subchanges[0] is canonical


@pytest.mark.parametrize('accepted', [True, False])
def test_run_moves_single(default_unidirectional_tis, accepted):
    # check that a single accepted step gives an MCStep with results that
    # match the expected active for accepted/rejected steps
    scheme = default_unidirectional_tis.scheme
    make_tis_trajectory = default_unidirectional_tis.make_tis_trajectory
    traj = {True: make_tis_trajectory(5),
            False: make_tis_trajectory(10)}[accepted]
    init_conds = scheme.initial_conditions_from_trajectories(traj)
    ensemble = scheme.network.sampling_ensembles[0]
    move = MockPathReversal(scheme, ensemble=ensemble)
    steplist = list(run_moves(init_conds, [move]))

    assert len(steplist) == 1
    step = steplist[0]
    assert isinstance(step, paths.MCStep)
    assert isinstance(step.change.canonical.mover, paths.PathReversalMover)
    assert step.change.accepted is accepted
    if accepted:
        assert step.active[ensemble].trajectory == traj.reversed
    else:
        assert step.active[ensemble] is init_conds[ensemble]


def test_run_moves_multiple(default_unidirectional_tis):
    scheme = default_unidirectional_tis.scheme
    traj = default_unidirectional_tis.make_tis_trajectory(10)
    make_trajectory = default_unidirectional_tis.make_trajectory
    init_conds = scheme.initial_conditions_from_trajectories(traj)
    ensemble = scheme.network.sampling_ensembles[2]

    moves = [
        # first move is force-accepted (always accepted anyway)
        MockForwardShooting(
            shooting_index=8,
            partial_traj=make_trajectory(8, 10),
            accepted=True,
            scheme=scheme,
            ensemble=ensemble
        ),
        # second move is rejected (bad ensemble)
        MockBackwardShooting(
            shooting_index=4,
            partial_traj=make_trajectory(3, 10),
            scheme=scheme,
            ensemble=ensemble
        ),
        # third move is force accepted
        MockBackwardShooting(
            shooting_index=9,
            partial_traj=make_trajectory(-1, 8).reversed,
            accepted=True,
            scheme=scheme,
            ensemble=ensemble
        )
    ]

    steps = list(run_moves(init_conds, moves))
    initial_traj = init_conds[ensemble].trajectory
    final_traj = steps[-1].active[ensemble].trajectory

    assert len(steps) == 3
    assert [step.change.accepted for step in steps] == [True, False, True]
    assert not initial_traj.is_correlated(final_traj)
