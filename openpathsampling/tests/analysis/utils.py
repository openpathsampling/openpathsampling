"""
Tools for generating sequences of steps for analysis tests.

The idea here is to have client code specify things like which snapshot
index is the shooting point and what partial trajectory from one-way
shooting should be, then we patch those things and use the internals of the
mover to generate the real data.

In this way, we use as much of the actual OPS machinery as possible,
ensuring, for example, that the details we return are correct.
"""


import numpy as np
import random
from unittest.mock import Mock, patch
from openpathsampling.tests.test_helpers import make_1d_traj

class AnalysisTestSetupError(Exception):
    """Raised for when an internal error occurs during test setup.

    These usually indicate a problem with test suite, not with the code
    itself.
    """
    pass

# TODO: add lower_bound support
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
        return random.choice(movers)

    try:
        signature = tuple(ensembles)
    except TypeError:
        signature = tuple([ensembles])
    for m in movers: print(m.ensemble_signature[0])
    sel = [m for m in movers if m.ensemble_signature[0] == signature]
    if len(sel) != 1:
        raise AnalysisTestSetupError(
            "expected 1 mover matching signature %s; found %d" %
            (signature, len(sel))
        )
    return sel[0]


class MockMove(object):
    def __init__(self, scheme, ensembles=None, group_name=None):
        self.scheme = scheme
        self.ensembles = ensembles
        self.group_name = group_name

    def mock_mover(self, mover):
        return mover

    def _generate_step_acceptance(self, inputs, accepted):
        for mover in random.shuffle(list(self.scheme[self.group_name])):
            change = self._do_move(mover, inputs)
            if change.accepted is accepted:
                return change

        looking_for = 'accepted' if accepted else 'rejected'
        raise AnalysisTestSetupError("unable to generate %s step" %
                                     looking_for)

    def accepted_step(self, inputs):
        return _generate_step_acceptance(inputs, accepted=True)

    def rejected_step(self, inputs):
        return _generate_step_acceptance(inputs, accepted=False)

    def _do_move(self, mover, inputs):
        mover = self.mock_mover(mover)
        change = mover.move(inputs)
        return change

    def __call__(self, inputs):
        mover = _select_by_input_ensembles(
            movers=self.scheme.movers[self.group_name],
            ensembles=self.ensembles
        )
        return self._do_move(mover, inputs)

class _MockSingleEnsembleMove(MockMove):
    def __init__(self, scheme, ensemble, group_name):
        super(_MockSingleEnsembleMove, self).__init__(
            scheme=scheme,
            ensembles=ensemble,
            group_name=group_name
        )

    @property
    def ensemble(self):
        return self.ensembles

class _MockOneWayShooting(_MockSingleEnsembleMove):
    def __init__(self, shooting_index, partial_traj, direction, scheme,
                 ensemble, group_name):
        super(_MockOneWayShooting, self).__init__(scheme=scheme,
                                                  ensemble=ensemble,
                                                  group_name=group_name)
        self.shooting_index = shooting_index
        self.partial_traj = partial_traj
        self.direction = direction

    def mock_move(self, mover):
        ...  # here's all the logic


class MockForwardShooting(_MockOneWayShooting):
    def __init__(self, shooting_index, partial_traj, scheme, ensemble=None,
                 group_name='shooting'):
        super(MockForwardShooting, self).__init__(
            shooting_index=shooting_index,
            partial_traj=partial_traj,
            direction='forward',
            scheme=scheme,
            ensemble=ensemble,
            group_name=group_name
        )

class MockBackwardShooting(_MockOneWayShooting):
    def __init__(self, shooting_index, partial_traj, scheme, ensemble=None,
                 group_name='shooting'):
        super(MockBackwardShooting, self).__init__(
            shooting_index=shooting_index,
            partial_traj=partial_traj,
            direction='backward',
            scheme=scheme,
            ensemble=ensemble,
            group_name=group_name
        )


class MockRepex(MockMove):
    def __init__(self, scheme, ensembles=None, group_name='repex'):
        super(MockRepex, self).__init__(scheme=scheme,
                                        ensembles=ensembles,
                                        group_name=group_name)

class MockPathReversal(_MockSingleEnsembleMove):
    def __init__(self, scheme, ensembles=None, group_name='pathreversal'):
        super(MockPathReversal, self).__init__(scheme=scheme,
                                               ensemble=ensembles,
                                               group_name=group_name)


def _do_single_step(init_conds, move, org_by_group):
    change = move(init_conds)
    if org_by_group:
        ...  # TODO: add this
    return change

def steps(init_conds, moves, org_by_group=True):
    for stepnum, move in enumerate(moves):
        change = _do_single_step(init_conds, move,
                                 org_by_group=org_by_group)
        new_conds = init_conds.apply_samples(change.results)
        step = paths.MCStep(
            simulation=None,
            mccycle=stepnum,
            active=new_conds,
            change=change
        )
        yield step
        init_conds = new_conds
