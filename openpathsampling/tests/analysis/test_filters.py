from openpathsampling.analysis.filters import *
import openpathsampling as paths
import numpy as np

try:
    from unittest.mock import Mock, patch
except ImportError:
    from mock import Mock, patch

from openpathsampling.tests.analysis.utils.mock_movers import (
    MockForwardShooting, MockBackwardShooting, MockRepex, MockPathReversal,
    MockTwoWayShooting, run_moves,
)

import pytest


### FIXTURES ###############################################################

@pytest.fixture
def scheme(default_unidirectional_tis):
    return default_unidirectional_tis.scheme


@pytest.fixture
def seven_steps_no_dynamics(default_unidirectional_tis):
    scheme = default_unidirectional_tis.scheme
    ens0, ens1, ens2 = scheme.network.sampling_ensembles
    t1 = default_unidirectional_tis.make_tis_trajectory(5)
    t2 = default_unidirectional_tis.make_tis_trajectory(10)
    init_conds = scheme.initial_conditions_from_trajectories([t1, t2])
    moves = [
        MockPathReversal(scheme, ensemble=ens0),  # ACC
        MockRepex(scheme, ensembles=[ens0, ens1]),  # ACC
        MockPathReversal(scheme, ensemble=ens2),  # REJ
        MockRepex(scheme, ensembles=[ens0, ens1]),  # ACC
        MockPathReversal(scheme, ensemble=ens0),  # ACC
        MockRepex(scheme, ensembles=[ens1, ens2]),  # REJ
        MockRepex(scheme, ensembles=[ens1, ens2]),  # REJ
    ]
    return run_moves(init_conds, moves)


@pytest.fixture
def seven_steps_trial_samples(seven_steps_no_dynamics):
    return sum([step.change.canonical.trials
                for step in seven_steps_no_dynamics], [])


@pytest.fixture
def shooting_and_pathrev_steps(default_unidirectional_tis):
    # 3 steps shooting (2 acc, 1 rej) and 1 step repex
    scheme = default_unidirectional_tis.scheme
    ens0, ens1, ens2 = scheme.network.sampling_ensembles
    t0 = default_unidirectional_tis.make_tis_trajectory(5)
    t1 = default_unidirectional_tis.make_tis_trajectory(10)
    init_conds = scheme.initial_conditions_from_trajectories([t0, t1])
    # t2: rejected trial in ens1 from shooting point ...
    p2 = ...
    # t3: accepted trial in ens1 from shooting point ...
    t3 = ...
    # t4: two-way shooting trial from shooting point ...
    t4 = ...
    # trials: (format: {replica: (traj, ens)})
    # step0 (initial):     {0: (t0, ens0), 1: (t0, ens1), 2: (t1, ens2)}
    # step1 (pathrev-acc): {0: (t0*, ens0), 1: (t0, ens1), 2: (t1, ens1)}
    # step2 (pathrev-acc): {0: (t0, ens0), 1: (t0, ens1), 2: (t1, ens1)}
    # step3 (shoot-rej):   {0: (t0, ens0), 1: (t2, ens1), 2: (t1, ens2)}
    # step4 (shoot-acc):   {0: (t0, ens0), 1: (t3, ens1), 2: (t1, ens2)}
    # step5 (twoway-acc):  {0: (t0, ens0), 1: (t4, ens1), 2: (t1, ens2)}
    # What these steps need to test:
    # 1. accepted and rejected shooting moves with same shooting point
    # 2. modified shooting point that differs from original
    # 3. includes non-trivial canonical mover (one-way shooting works)
    # 4. includes non-shooting moves as well
    moves = [
        MockPathReversal(scheme, ensemble=ens0),  # ACC
        MockPathReversal(scheme, ensemble=ens0),  # ACC
        MockForwardShooting(shooting_index=4,
                            partial_traj=p2,
                            scheme=scheme,
                            ensemble=ens1),
        MockForwardShooting(shooting_index=4,
                            partial_traj=p3,
                            scheme=scheme,
                            ensemble=ens1),
        MockTwoWayShooting(shooting_index=4,
                           new_traj=t4,
                           modified_snapshot=...,
                           direction="forward-first",
                           scheme=scheme,
                           ensemble=ens1,
                           accepted=True),
    ]
    return run_moves(init_conds, moves)



### TESTS ##################################################################

##### generic filter logic #################################################

class StringFilter(GenericFilter):
    FILTER_TYPE = str


@pytest.fixture
def strings():
    return ['foo', 'bar', 'baz', 'qux', 'quux']


def test_filter_condition(strings):
    has_a = StringFilter(condition=lambda s: 'a' in s, name='has_a')
    assert list(has_a(strings)) == ['bar', 'baz']


def test_or_filter_condition(strings):
    has_a = StringFilter(condition=lambda s: 'a' in s, name='has_a')
    len_4 = StringFilter(condition=lambda s: len(s) == 4, name='len_4')
    filt = has_a | len_4
    assert list(filt(strings)) == ['bar', 'baz', 'quux']


def test_and_filter_condition(strings):
    has_q = StringFilter(condition=lambda s: 'q' in s, name='has_q')
    len_4 = StringFilter(condition=lambda s: len(s) == 4, name='len_4')
    filt = has_q & len_4
    assert list(filt(strings)) == ['quux']


def test_invert_filter_condition(strings):
    has_q = StringFilter(condition=lambda s: 'q' in s, name='has_q')
    filt = ~has_q
    assert list(filt(strings)) == ['foo', 'bar', 'baz']


def test_sub_filter_condition(strings):
    len_3 = StringFilter(condition=lambda s: len(s) == 3, name='len_3')
    has_b = StringFilter(condition=lambda s: 'b' in s, name='has_b')
    filt = len_3 - has_b
    assert list(filt(strings)) == ['foo', 'qux']


def test_xor_filter_condition(strings):
    len_3 = StringFilter(condition=lambda s: len(s) == 3, name='len_3')
    has_q = StringFilter(condition=lambda s: 'q' in s, name='has_q')
    filt = len_3 ^ has_q
    assert list(filt(strings)) == ['foo', 'bar', 'baz', 'quux']


@pytest.mark.parametrize('returns_iterable', [True, False])
def test_extractor(strings, returns_iterable):
    extractor_func = {
        True: lambda x: x[:2],
        False: lambda x: x[0]
    }[returns_iterable]
    expected = {
        False: ['f', 'b', 'b', 'q', 'q'],
        True: ['fo', 'ba', 'ba', 'qu', 'qu'],
    }[returns_iterable]

    extractor = Extractor(extractor_func, name='string_ext',
                          returns_iterable=returns_iterable)
    extracted = [extractor(s) for s in strings]
    assert extracted == expected


def test_extractor_using(strings):
    pytest.skip()


##### step filters #########################################################

@pytest.mark.parametrize('input_type', ['mover', 'class_name', 'class'])
def test_canonical_mover_step_filter(input_type, seven_steps_no_dynamics,
                                     scheme):
    input_param = {
        'mover': scheme.movers['pathreversal'][0],
        'class_name': 'PathReversalMover',
        'class': paths.PathReversalMover
    }[input_type]
    expected = {
        'mover': [0, 4],
        'class_name': [0, 2, 4],
        'class': [0, 2, 4]
    }[input_type]

    filt = canonical_mover(input_param)
    filtered_steps = list(filt(seven_steps_no_dynamics))
    assert len(filtered_steps) == len(expected)
    assert [s.mccycle for s in filtered_steps] == expected


def test_canonical_mover_step_filter_input_step(seven_steps_no_dynamics):
    step = next(seven_steps_no_dynamics)
    with pytest.raises(TypeError, match="mean to use the extractor"):
        _ = canonical_mover(step)


def test_canonical_mover_step_filter_bad_input():
    with pytest.raises(TypeError, match="does not appear to be"):
        _ = canonical_mover(5)


def test_trial_replica_step_filter(seven_steps_no_dynamics):
    filt = trial_replica(0)
    filtered_steps = list(filt(seven_steps_no_dynamics))
    assert len(filtered_steps) == 4
    assert [s.mccycle for s in filtered_steps] == [0, 1, 3, 4]


def test_trial_ensemble_step_filter(scheme, seven_steps_no_dynamics):
    ens0 = scheme.network.sampling_ensembles[0]
    filt = trial_ensemble(ens0)
    filtered_steps = list(filt(seven_steps_no_dynamics))
    assert len(filtered_steps) == 4
    assert [s.mccycle for s in filtered_steps] == [0, 1, 3, 4]


def test_rejected_steps_filter(seven_steps_no_dynamics):
    filtered_steps = list(rejected_steps(seven_steps_no_dynamics))
    assert len(filtered_steps) == 3
    assert [s.mccycle for s in filtered_steps] == [2, 5, 6]


def test_accepted_steps_filter(seven_steps_no_dynamics):
    filtered_steps = list(accepted_steps(seven_steps_no_dynamics))
    assert len(filtered_steps) == 4
    assert [s.mccycle for s in filtered_steps] == [0, 1, 3, 4]


def test_all_steps_filter(seven_steps_no_dynamics):
    filtered_steps = list(all_steps(seven_steps_no_dynamics))
    assert len(filtered_steps) == 7
    assert [s.mccycle for s in filtered_steps] == list(range(7))


##### sample filters #######################################################

def test_all_samples_filter(seven_steps_trial_samples):
    samples = list(all_samples(seven_steps_trial_samples))
    assert len(samples) == 11


def test_ensemble_sample_filter(scheme, seven_steps_trial_samples):
    ens0, ens1, ens2 = scheme.network.sampling_ensembles
    expected = {ens0: 4, ens1: 4, ens2: 3}
    for ens in [ens0, ens1, ens2]:
        ens_filter = ensemble(ens)
        filtered = list(ens_filter(seven_steps_trial_samples))
        assert len(filtered) == expected[ens]


def test_replica_sample_filter(seven_steps_trial_samples):
    expected = {0: 4, 1: 4, 2: 3}
    for rep in [0, 1, 2]:
        rep_filter = replica(rep)
        filtered = list(rep_filter(seven_steps_trial_samples))
        assert len(filtered) == expected[rep]


def test_minus_ensemble_sample_filter():
    pytest.skip()


def test_ms_outer_ensemble_sample_filter():
    pytest.skip()


def test_sampling_ensemble_sample_filter():
    pytest.skip()


##### extractors ###########################################################

def _map_replicas_and_ensembles(sample_set):
    return {s.replica: s.ensemble for s in sample_set}


def test_active_samples_extractor(scheme, seven_steps_no_dynamics):
    ens0, ens1, ens2 = scheme.network.sampling_ensembles
    all_expected = [
        {0: ens0, 1: ens1, 2: ens2},
        {0: ens1, 1: ens0, 2: ens2},
        {0: ens1, 1: ens0, 2: ens2},
        {0: ens0, 1: ens1, 2: ens2},
        {0: ens0, 1: ens1, 2: ens2},
        {0: ens0, 1: ens1, 2: ens2},
        {0: ens0, 1: ens1, 2: ens2},
    ]
    for step, expected in zip(seven_steps_no_dynamics, all_expected):
        found = _map_replicas_and_ensembles(step.active)
        assert found == expected


def test_trial_samples_extractor(scheme, seven_steps_no_dynamics):
    ens0, ens1, ens2 = scheme.network.sampling_ensembles
    all_expected = [
        {0: ens0},
        {0: ens1, 1: ens0},
        {2: ens2},
        {0: ens0, 1: ens1},
        {0: ens0},
        {1: ens2, 2: ens1},
        {1: ens2, 2: ens1},
    ]
    for step, expected in zip(seven_steps_no_dynamics, all_expected):
        found = _map_replicas_and_ensembles(step.change.trials)
        assert found == expected


def test_shooting_steps_filter(shooting_and_pathrev_steps):
    shooting = list(shooting_steps(shooting_and_pathrev_steps))
    assert len(shooting) == 3
    pytest.skip()
    for step in shooting:
        assert isinstance(step.change.canonical.mover,
                          (paths.ForwardFirstTwoWayShooting,
                           paths.BackwardFirstTwoWayShooing,
                           paths.ForwardShooting,
                           paths.BackwardShooting))


def test_shooting_points_extractor(shooting_and_pathrev_steps):
    pytest.skip()


def test_modified_shooting_points_extractor(shooting_and_pathrev_steps):
    pytest.skip()


def test_canonical_movers_extractor(scheme, seven_steps_no_dynamics):
    pytest.skip()
