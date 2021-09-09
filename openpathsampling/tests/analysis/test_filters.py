from openpathsampling.analysis.filters import *
import openpathsampling as paths
import numpy as np

from unittest import Mock, patch

import pytest

### TESTS ##################################################################

##### step filters #########################################################

@pytest.mark.parametrize('input_type', ['mover', 'group_name', 'class'])
def test_canonical_mover_step_filter(input_type):
    pytest.skip()

def test_trial_replica_step_filter():
    pytest.skip()

def test_trial_ensemble_step_filter():
    pytest.skip()

def test_rejected_steps_filter():
    pytest.skip()

def test_accepted_steps_filter():
    pytest.skip()

def test_all_steps_filter():
    pytest.skip()

def test_all_samples_filter():
    pytest.skip()

def test_ensemble_sample_filter():
    pass

def test_replica_sample_filter():
    pass

def test_minus_ensemble_sample_filter():
    pass

def test_ms_outer_ensemble_sample_filter():
    pass

def test_sampling_ensemble_sample_filter():
    pass

def test_active_samples_extractor():
    pass

def test_trial_samples_extractor():
    pass

def test_shooting_steps_filter():
    pass

def test_shooting_points_extractor():
    pass

def test_modified_shooting_points_extractor():
    pass

def test_canonical_movers_extractor():
    pass


