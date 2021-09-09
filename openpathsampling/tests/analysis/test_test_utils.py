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
    print([m.ensemble_signature for m in movers])
    if ensemble == 'ensemble':
        ensemble = frozenset([movers[0]])

    selected = _select_by_input_ensembles(movers, ensemble)
    assert selected in movers  # this is enough for ensemble=None
    pytest.skip()

def test_wrap_one_way_shooting_mover(scheme):
    pytest.skip()

def test_wrap_repex_move(scheme):
    pytest.skip()
