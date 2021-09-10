"""
Tools for generating sequences of steps for analysis tests.

The idea here is to have client code specify things like which snapshot
index is the shooting point and what partial trajectory from one-way
shooting should be, then we patch those things and use the internals of the
mover to generate the real data.

In this way, we use as much of the actual OPS machinery as possible,
ensuring, for example, that the details we return are correct.

Note that for the classes here, each instance represents a single move
(single MC step). This is unlike the PathMover objects in OPS, which
represents an individual type of move that can be reused for many steps.
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
def make_trajectory(cv_max, lower_bound=None):
    """Create a trajectory with maximum x at cv_max < x < cv_max + 0.1
    """
    do_decreasing = False
    if lower_bound is None:
        lower_bound = -0.1
        do_decreasing = True
    increasing = list(np.arange(lower_bound, cv_max + 0.01, 0.1))
    if do_decreasing and cv_max < 1.0:
        decreasing = list(reversed(increasing))[1:]
    else:
        decreasing = []
    array = np.array(increasing + decreasing) + np.random.random() * 0.1
    return make_1d_traj(array)

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
            "expected 1 mover matching signature %s; found %d. Allowed: %s"
            % (signature, len(sel), [m.ensemble_signature[0] for m in movers])
        )
    return sel[0]

def _run_patched(mover, patches, inputs):
    for patch in patches:
        patch.start()

    change = mover.move(inputs)

    for patch in patches:
        patch.stop()
    return change


class MockRandomChoiceMover(object):
    def __init__(self, random_mover, change):
        self.random_mover = random_mover
        self.change = change

    def patches(self):
        mover = self.change.mover
        try:
            idx = self.random_mover.submovers.index(mover)
        except ValueError:
            raise AnalysisTestSetupError(
                "mover %s not found as submover of mover %s" %
                (self.random_mover, mover)
            )
        rng_mock = Mock(choice=Mock(return_value=idx))
        patches = [
            patch.object(self.random_mover, '_rng', rng_mock),
            patch.object(mover, 'move', Mock(return_value=self.change))
        ]
        return patches

    def __call__(self, inputs):
        patches = self.patches()
        return _run_patched(self.random_mover, patches, inputs)


class MockMove(object):
    def __init__(self, scheme, ensembles=None, group_name=None):
        self.scheme = scheme
        self.ensembles = ensembles
        self.group_name = group_name

    def wrap_org_by_group(self, change, inputs):
        # extract the root_mover (selects type of move) and the
        # group_selector (selects a specific move within the move type)
        root_mover = self.scheme.root_mover
        group_selectors = root_mover.submovers
        group_selector = [g for g in group_selectors
                          if change.mover in g.submovers]
        if len(group_selector) != 1:
            raise AnalysisTestSetupError(
                "expected 1 group containing the mover %s; found %d" %
                (change.mover, len(group_selector))
            )
        group_selector = group_selector[0]

        # make a move change for the inner step (selecting which mover
        # within the move type)
        inner_mock = MockRandomChoiceMover(group_selector, change)
        inner_change = inner_mock(inputs)

        # make a move change for the outer step (selecting which move type
        # to do from the root_mover)
        outer_mock = MockRandomChoiceMover(root_mover, inner_change)
        outer_change=  outer_mock(inputs)
        return outer_change

    def patches(self, mover):
        return []

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
        patches = self.patches(mover)
        return _run_patched(mover, patches, inputs)

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
                 ensemble, accepted, group_name):
        super(_MockOneWayShooting, self).__init__(scheme=scheme,
                                                  ensemble=ensemble,
                                                  group_name=group_name)
        self.shooting_index = shooting_index
        self.partial_traj = partial_traj
        self.direction = direction
        self.accepted = accepted

    def patches(self, mover):
        selector_mock = Mock(pick=Mock(return_value=self.shooting_index))
        if self.accepted is not None:
            value = 1.0 if accepted else 0.0
            selector_mock.probability_ratio = Mock(return_value=value)

        selector_patch = patch.object(mover, 'selector', selector_mock)

        if self.direction == 'forward':
            traj_patch = patch.object(mover, '_make_forward_trajectory',
                                      Mock(return_value=self.partial_traj))
        elif self.direction == 'backward':
            traj_patch = patch.object(mover, '_make_backward_trajectory',
                                      Mock(return_value=self.partial_traj))
        else:
            raise AnalysisTestSetupError("Invalid direction: %s" %
                                         self.direction)

        return [selector_patch, traj_patch]


class MockForwardShooting(_MockOneWayShooting):
    def __init__(self, shooting_index, partial_traj, scheme, ensemble=None,
                 accepted=None, group_name='shooting'):
        super(MockForwardShooting, self).__init__(
            shooting_index=shooting_index,
            partial_traj=partial_traj,
            direction='forward',
            scheme=scheme,
            ensemble=ensemble,
            accepted=accepted,
            group_name=group_name
        )


class MockBackwardShooting(_MockOneWayShooting):
    def __init__(self, shooting_index, partial_traj, scheme, ensemble=None,
                 accepted=None, group_name='shooting'):
        super(MockBackwardShooting, self).__init__(
            shooting_index=shooting_index,
            partial_traj=partial_traj,
            direction='backward',
            scheme=scheme,
            ensemble=ensemble,
            accepted=accepted,
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
        change = move.wrap_org_by_group(change, init_conds)
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
