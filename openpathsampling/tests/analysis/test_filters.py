from openpathsampling.analysis.filters import *
import openpathsampling as paths
import numpy as np

from unittest.mock import Mock, patch

from openpathsampling.tests.analysis.utils import (
    MockForwardShooting, MockBackwardShooting, MockRepex, MockPathReversal,
    make_tis_trajectory, make_trajectory, run_moves
)

import pytest


### FIXTURES ###############################################################

@pytest.fixture
def seven_steps_no_dynamics(scheme):
    ens0, ens1, ens2 = scheme.network.sampling_ensembles
    t1 = make_tis_trajectory(0.5)
    t2 = make_tis_trajectory(1.0)
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


def test_extractor(strings):
    pytest.skip()

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
    pytest.skip()

def test_ensemble_sample_filter(scheme, seven_steps_trial_samples):
    pytest.skip()

def test_replica_sample_filter(seven_steps_trial_samples):
    pytest.skip()

def test_minus_ensemble_sample_filter():
    pytest.skip()

def test_ms_outer_ensemble_sample_filter():
    pytest.skip()

def test_sampling_ensemble_sample_filter():
    pytest.skip()

##### extractors ###########################################################

def test_active_samples_extractor():
    pytest.skip()

def test_trial_samples_extractor():
    pytest.skip()

def test_shooting_steps_filter():
    pytest.skip()

def test_shooting_points_extractor():
    pytest.skip()

def test_modified_shooting_points_extractor():
    pytest.skip()

def test_canonical_movers_extractor():
    pytest.skip()
