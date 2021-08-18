from openpathsampling.analysis.filters import *
import openpathsampling as paths

import pytest

@pytest.fixture
def tis_network():
    cv = paths.FunctionCV("x", lambda x: x)
    state_A = paths.CVDefinedVolume(cv, float("-inf"), 0)
    state_B = paths.CVDefinedVolume(cv, 1.0, float("inf"))
    interfaces = paths.VolumeInterfaceSet(cv, float("-inf"), [0.0, 0.1, 0.2])
    network = paths.MISTISNetwork([(state_A, interfaces, state_B)])
    return network

@pytest.fixture
def scheme(tis_network):
    engine = paths.engines.NoEngine({'n_spatial': 3, 'n_atoms': 1})
    scheme = paths.DefaultScheme(network=tis_network, engine=engine)
    return scheme

### HELPERS ################################################################

def _select_by_input_signature(movers, signature):
    sel = [m for m in movers if m.ensemble_signature == signature]
    assert len(sel) == 1
    return sel[0]

def make_shooting_move(scheme, accepted=None, ensemble=None):
    pass

def make_repex_move(scheme, accepted=None, ensembles=None):
    pass

def make_msouter_move(scheme, accepted=None):
    pass

def make_minus_move(scheme, accepted=None):
    pass

def make_steps(moves):
    pass


### TESTS ##################################################################

def test_canonical_mover_step_filter():
    pass

def test_trial_replica_step_filter():
    pass

def test_trial_ensemble_step_filter():
    pass

def test_rejected_steps_filter():
    pass

def test_accepted_steps_filter():
    pass

def test_all_steps_filter():
    pass

def test_all_samples_filter():
    pass

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


