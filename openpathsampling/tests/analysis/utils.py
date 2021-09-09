"""
Tools for generating sequences of steps for analysis tests.


"""
import numpy as np
from openpathsampling.tests.test_helpers import make_1d_traj

def make_trajectory(cv_max):
    """Create a trajectory with maximum x at cv_max < x < cv_max + 0.1
    """
    increasing = list(range(-1, int(cv_max * 10) + 1))
    if cv_max < 1.0:
        decreasing = list(reversed(increasing))[1:]
    else:
        decreasing = []
    array = np.array(increasing + decreasing) + np.random.random()
    return make_1d_traj(array * 0.1)

def _select_by_input_ensembles(movers, ensembles):
    # quick return if we don't actually care about which mover we use
    if ensembles is None:
        return np.random.choice(movers)

    try:
        signature = frozenset(ensembles)
    except TypeError:
        signature = frozenset([ensembles])
        sel = [m for m in movers if m.ensemble_signature == signature]
        if len(sel) != 1:
            raise RuntimeError("Error is test setup: expected 1 mover "
                               "matching signature; found %d" % len(sel))
        return sel[0]

def wrap_one_way_shooting_move(scheme, ensemble=None, shooting_index=None):
    pass

def wrap_repex_move(scheme, accepted=None, ensembles=None):
    if ensembles is not None:
        mover = _select_by_input_signature(
            movers=scheme.groups['repex'],
            signature=frozenset(ensembles)
        )
    else:
        mover = np.random.choice(scheme.groups['repex'])

    pass

def wrap_msouter_move(scheme, accepted=None):
    pass

def wrap_minus_move(scheme, accepted=None):
    pass

def make_steps(moves):
    pass



