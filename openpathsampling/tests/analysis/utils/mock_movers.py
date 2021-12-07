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

import random
try:
    from unittest.mock import Mock, patch
except ImportError:  # pragma: no cover
    from mock import Mock, patch

import openpathsampling as paths


def _get_only(iterable, condition, error_msg):
    possibilities = [item for item in iterable if condition(item)]
    if len(possibilities) != 1:
        msg = "expected 1 {error_msg}; found {nfound}: {found}".format(
            error_msg=error_msg,
            nfound=len(possibilities),
            found=possibilities
        )
        raise AnalysisTestSetupError(msg)
    return possibilities[0]


class AnalysisTestSetupError(Exception):
    """Raised for when an internal error occurs during test setup.

    These usually indicate a problem with test suite, not with the code
    itself.
    """


def _select_by_input_ensembles(movers, ensembles):
    """Select a mover by input ensembles.

    Pick a random mover if ensembles is ``None``.
    """
    # quick return if we don't actually care about which mover we use
    if ensembles is None:
        return random.choice(movers)

    try:
        signature = tuple(ensembles)
    except TypeError:
        signature = tuple([ensembles])

    mover = _get_only(
        iterable=movers,
        condition=lambda m: set(m.ensemble_signature[0]) == set(signature),
        error_msg="mover matching signature{sig}".format(sig=signature)
    )
    return mover


def _run_patched(mover, patches, inputs):
    """Run the move with the given patches for the given inputs"""
    for p in patches:
        p.start()

    change = mover.move(inputs)

    for p in patches:
        p.stop()
    return change


class MockRandomChoiceMover(object):
    """Mock random choice mover to give the desired subchange.

    This differs from the regular mover mocks in that it is only created
    when we have already determined the subchange that will be generated.

    Parameters
    ----------
    random_mover: class:`.RandomChoiceMover`
        the mover to mock the behavior or
    change : :class:`.MoveChange`
        the change that the submover would return
    """
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
    """Abstract object for mocking a single step of a path simulation.
    """
    def __init__(self, scheme, ensembles=None, group_name=None):
        self.scheme = scheme
        self.ensembles = ensembles
        self.group_name = group_name

    def wrap_org_by_group(self, change, inputs):
        """Wraps as if using an OrganizeByMoveGroup strategy
        """
        # extract the root_mover (selects type of move) and the
        # group_selector (selects a specific move within the move type)
        root_mover = self.scheme.root_mover
        group_selector = _get_only(
            iterable=root_mover.submovers,
            condition=lambda g: change.mover in g.submovers,
            error_msg="group containing the mover {mover}".format(
                mover=change.mover
            )
        )

        # make a move change for the inner step (selecting which mover
        # within the move type)
        inner_mock = MockRandomChoiceMover(group_selector, change)
        inner_change = inner_mock(inputs)

        # make a move change for the outer step (selecting which move type
        # to do from the root_mover)
        outer_mock = MockRandomChoiceMover(root_mover, inner_change)
        outer_change = outer_mock(inputs)
        return outer_change

    def patches(self, mover):
        """Patches to be used by this mover

        Subclasses will override this when they need to control the mover's
        behavior.

        Parameters
        ----------
        mover : :class:`.PathMover`
            the path mover to patch

        Returns
        -------
        List[unittest.patch]
            patch objects to use during the move
        """
        return []

    def safety_checks(self, change, inputs):
        """Safety checks to ensure reasonable output.

        Override in subclasses.

        Parameters
        ----------
        change : :class:`.MoveChange`
            move change from this mocked step
        inputs : :class:`.SampleSet`
            original input to this mover
        """
        pass

    def _do_move(self, mover, inputs):
        patches = self.patches(mover)
        return _run_patched(mover, patches, inputs)

    def __call__(self, inputs):
        mover = _select_by_input_ensembles(
            movers=self.scheme.movers[self.group_name],
            ensembles=self.ensembles
        )
        change = self._do_move(mover, inputs)
        self.safety_checks(change, inputs)
        return change


class _MockSingleEnsembleMove(MockMove):
    def __init__(self, scheme, ensemble, group_name):
        if ensemble is not None:
            ensemble = [ensemble]
        super(_MockSingleEnsembleMove, self).__init__(
            scheme=scheme,
            ensembles=ensemble,
            group_name=group_name
        )

    @property
    def ensemble(self):
        return self.ensembles[0]


class _MockOneWayShooting(_MockSingleEnsembleMove):
    """Abstract mock for one-way shooting.

    Core logic implemented in here (see methods ``patches`` and
    ``safety_checks``). For usage information, see MockForwardShooting and
    MockBackwardShooting.
    """
    def __init__(self, shooting_index, partial_traj, direction, scheme,
                 ensemble, accepted, group_name):
        super(_MockOneWayShooting, self).__init__(scheme=scheme,
                                                  ensemble=ensemble,
                                                  group_name=group_name)
        self.shooting_index = shooting_index
        self.partial_traj = partial_traj
        self.direction = direction
        self.accepted = accepted

    def _generate_mock(self, snapshot, running):
        return paths.Trajectory([snapshot] + self.partial_traj)

    def patches(self, mover):
        if self.direction == 'forward':
            shoot_type = paths.ForwardShootMover
        elif self.direction == 'backward':
            shoot_type = paths.BackwardShootMover
        else:  # pragma: no cover
            raise AnalysisTestSetupError("Invalid direction: %s" %
                                         self.direction)

        shooter = _get_only(
            iterable=mover.submovers,
            condition=lambda m: isinstance(m, shoot_type),
            error_msg=(
                "submover of type {shoot_type};".format(
                    shoot_type=shoot_type.__class__.__name__)
                )
            )

        # patch for the choice of fwd vs bkwd shooter
        direction_idx = mover.submovers.index(shooter)
        rng_mock = Mock(choice=Mock(return_value=direction_idx))
        direction_patch = patch.object(mover, '_rng', rng_mock)

        # patch for the partial trajectory returned by dynamics
        traj_patch = patch.object(shooter.engine, 'generate',
                                  self._generate_mock)

        # patch for the shooting point selector
        pick_patch = patch.object(shooter.selector, 'pick',
                                  Mock(return_value=self.shooting_index))
        if self.accepted is not None:
            value = 1.0 if self.accepted else 0.0
            probability_ratio_patch = [patch.object(shooter.selector,
                                                    'probability_ratio',
                                                    return_value=value)]
        else:
            probability_ratio_patch = []

        return ([direction_patch, traj_patch, pick_patch]
                + probability_ratio_patch)

    def _check_shooting_point_in_trial(self, change, inputs):
        # was the shooting point included in the trial trajectory? (we can't
        # know this until after we create the move change)
        shooting_snapshot = change.canonical.details.shooting_snapshot
        trial_traj = change.canonical.trials[0].trajectory
        if shooting_snapshot not in trial_traj:  # pragma: no cover
            # This should not be possible, because the mock mover adds the
            # shooting point. Safety against weird failures, though.
            raise AnalysisTestSetupError("shooting point not included in "
                                         "partial trajectory")

    def _check_forced_accepted_reasonable(self, change, inputs):
        trial = change.canonical.trials[0]
        if self.accepted and not trial.ensemble(trial.trajectory):
            if change.accepted:  # pragma: no cover
                raise AnalysisTestSetupError(
                    "Something when wrong in the mock mover. A step that "
                    "can not be accepted was accepted."
                )
            else:
                raise AnalysisTestSetupError(
                    "Step tried to force acceptance on a trial that can "
                    "not be accepted."
                )

    def safety_checks(self, change, inputs):
        self._check_shooting_point_in_trial(change, inputs)
        self._check_forced_accepted_reasonable(change, inputs)


class MockForwardShooting(_MockOneWayShooting):
    """Mock a forward shooting move.

    Parameters
    ----------
    shooting_index : int
        the index of the shooting point in the input trajectory
    partial_traj : :class:`.Trajectory`
        The trajectory to mock as coming from the engine. This is the
        trajectory *without* the shooting point -- new frames only.
    scheme : :class:`.MoveScheme`
        the move scheme used in the mock simulation
    ensemble : :class:`.Ensemble` or None
        the ensemble to shoot from. If ``None``, a random mover is selected
        from the movers in ``group_name``.
    accepted : bool or None
        if ``True`` or ``False``, force the acceptance to by true or false,
        respectively. If ``None`` (default) use the internal logic to
        determine acceptance. NOTE: This only bypasses the Monte Carlo
        proposal probability calculation. Moves that are lot legally
        possible to accept will raise an AnalysisTestSetupError if you try
        to force them to be accepted.
    group_name : str
        the name of the groups to select movers from (default
        ``'shooting'``).
    """
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
    """Mock a backward shooting move.

    Parameters
    ----------
    shooting_index : int
        the index of the shooting point in the input trajectory
    partial_traj : :class:`.Trajectory`
        The trajectory to mock as coming from the engine. This is the
        trajectory *without* the shooting point -- new frames only.
    scheme : :class:`.MoveScheme`
        the move scheme used in the mock simulation
    ensemble : :class:`.Ensemble` or None
        the ensemble to shoot from. If ``None``, a random mover is selected
        from the movers in ``group_name``.
    accepted : bool or None
        if ``True`` or ``False``, force the acceptance to by true or false,
        respectively. If ``None`` (default) use the internal logic to
        determine acceptance. NOTE: This only bypasses the Monte Carlo
        proposal probability calculation. Moves that are lot legally
        possible to accept will raise an AnalysisTestSetupError if you try
        to force them to be accepted.
    group_name : str
        the name of the groups to select movers from (default
        ``'shooting'``).
    """

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
    """Mock a replica exchange move.

    Parameters
    ----------
    scheme : :class:`.MoveScheme`
        the move scheme used in the mock simulation
    ensembles : List[:class:`.Ensemble`] or None
        The input ensembles for the replica exchange move. If ``None``, a
        random mover is selected.
    group_name : str
        the name of the groups to select movers from (default
        ``'repex'``).
    """
    def __init__(self, scheme, ensembles=None, group_name='repex'):
        super(MockRepex, self).__init__(scheme=scheme,
                                        ensembles=ensembles,
                                        group_name=group_name)


class MockPathReversal(_MockSingleEnsembleMove):
    """Mock a path reversal  move.

    Parameters
    ----------
    scheme : :class:`.MoveScheme`
        the move scheme used in the mock simulation
    ensemble : List[:class:`.Ensemble`] or None
        The input ensemble for the path reversal move. If ``None``, a
        random mover is selected.
    group_name : str
        the name of the groups to select movers from (default
        ``'pathreversal'``).
    """
    def __init__(self, scheme, ensemble=None, group_name='pathreversal'):
        super(MockPathReversal, self).__init__(scheme=scheme,
                                               ensemble=ensemble,
                                               group_name=group_name)


def _do_single_step(init_conds, move, org_by_group):
    change = move(init_conds)
    if org_by_group:
        change = move.wrap_org_by_group(change, init_conds)
    return change


def run_moves(init_conds, moves, org_by_group=True):
    """Run a sequence of moves.

    Once the individual moves have been mocked up, this runs them to iterate
    over the steps that you would get from the simulation.

    Since this is a generator, you may want to wrap it in a list when using
    it.

    Parameters
    ----------
    init_conds : :class:`.SampleSet`
        initial conditions for this run
    moves : List[:class:`.MockMove`]
        the moves to run, in the order they will be applied
    org_by_group : bool
        whether to wrap each resulting change to include the structural move
        changes that come from the :class:`.OrganizeByMoveGroupStrategy`.

    Yields
    ------
    :class:`.MCStep`
        the step for each mock move
    """
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
