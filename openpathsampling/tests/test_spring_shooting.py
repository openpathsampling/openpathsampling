import pytest
import openpathsampling as paths
from openpathsampling import MoveChange, Details
from openpathsampling.high_level.move_scheme import MoveScheme
from openpathsampling.tests.test_movestrategy import MoveStrategyTestSetup
from openpathsampling.pathmover import SampleNaNError, SampleMaxLengthError
from openpathsampling.tests.test_helpers import (make_1d_traj,
                                                 assert_items_equal,
                                                 CalvinistDynamics,
                                                 NaNEngine)

from openpathsampling.pathmovers.spring_shooting import (
    SpringShootingSelector, SpringMover,
    ForwardSpringMover, BackwardSpringMover,
    SpringShootingMover, SpringShootingStrategy,
    SpringShootingMoveScheme)


class FakeStep(object):
    def __init__(self, details):
        self.change = MoveChange(details=details)


class SelectorTest(object):
    def setup(self):
        self.mytraj = make_1d_traj(coordinates=[-0.5, 0.1, 0.2, 0.3, 0.5],
                                   velocities=[1.0, 1.0, 1.0, 1.0, 1.0])
        self.dyn = CalvinistDynamics([-0.5, -0.4, -0.3, -0.2, -0.1,
                                      0.1, 0.2, 0.3, 0.4, 0.5])
        self.dyn.initialized = True
        self.initial_guess = 3
        self.ens = paths.LengthEnsemble(5)
        self.ges = paths.SampleSet(paths.Sample(
            replica=0, trajectory=self.mytraj, ensemble=self.ens
        ))


class TestSpringShootingSelector(SelectorTest):

    @staticmethod
    def test_neg_delta():
        with pytest.raises(ValueError):
            SpringShootingSelector(delta_max=-1,
                                   k_spring=0.1)

    @staticmethod
    def test_0_delta():
        with pytest.raises(ValueError):
            SpringShootingSelector(delta_max=0,
                                   k_spring=0.1)

    @staticmethod
    def test_pos_delta():
        delta_max = 1
        sel = SpringShootingSelector(delta_max=delta_max,
                                     k_spring=0.1)
        assert sel.delta_max == delta_max
        assert len(sel._fw_prob_list) == 2*delta_max+1
        assert len(sel._bw_prob_list) == 2*delta_max+1

    @staticmethod
    def test_neg_k():
        with pytest.raises(ValueError):
            SpringShootingSelector(delta_max=1,
                                   k_spring=-1)

    @staticmethod
    def test_0_k():
        sel = SpringShootingSelector(delta_max=1,
                                     k_spring=0)
        ref_list = [1.0, 1.0, 1.0]
        assert sel.k_spring == 0
        assert sel._fw_prob_list == ref_list
        assert sel._bw_prob_list == ref_list

    @staticmethod
    def test_pos_k():
        sel = SpringShootingSelector(delta_max=1,
                                     k_spring=1)
        assert sel.k_spring == 1

    @staticmethod
    def test_sanity_breaking_fw():
        sel = SpringShootingSelector(delta_max=1,
                                     k_spring=1)
        sel._fw_prob_list = [0, 0, 0]
        with pytest.raises(RuntimeError):
            sel.check_sanity()

    @staticmethod
    def test_sanity_breaking_bw():
        sel = SpringShootingSelector(delta_max=1,
                                     k_spring=1)
        sel._bw_prob_list = [0, 0, 0]

        with pytest.raises(RuntimeError):
            sel.check_sanity()

    @staticmethod
    def test_sanity_breaking_total():
        sel = SpringShootingSelector(delta_max=1,
                                     k_spring=1)
        sel._total_bias = sum([0, 0, 0])
        with pytest.raises(RuntimeError):
            sel.check_sanity()

    @staticmethod
    def test_probability_ratio():
        sel = SpringShootingSelector(delta_max=1,
                                     k_spring=1)
        sel._acceptable_snapshot = True
        ratio = sel.probability_ratio(None, None, None)
        assert ratio == 1.0
        sel._acceptable_snapshot = False
        ratio = sel.probability_ratio(None, None, None)
        assert ratio == 0.0

    def test_sanity_breaking_fw_pick(self):
        sel = SpringShootingSelector(delta_max=1,
                                     k_spring=1)
        sel._fw_prob_list = [0, 0, 0]
        with pytest.raises(RuntimeError):
            sel.pick(trajectory=self.mytraj, direction='forward')

    def test_sanity_breaking_bw_pick(self):
        sel = SpringShootingSelector(delta_max=1,
                                     k_spring=1)
        sel._bw_prob_list = [0, 0, 0]
        with pytest.raises(RuntimeError):
            sel.pick(trajectory=self.mytraj, direction='forward')

    def test_sanity_breaking_total_pick(self):
        sel = SpringShootingSelector(delta_max=1,
                                     k_spring=1)
        sel._total_bias = sum([0, 0, 0])
        with pytest.raises(RuntimeError):
            sel.pick(trajectory=self.mytraj, direction='forward')

    @staticmethod
    def test_initial_guess():
        sel = SpringShootingSelector(delta_max=1, k_spring=1, initial_guess=12)
        assert sel.trial_snapshot == 12
        assert sel.previous_snapshot == 12

    def test_pick_direction(self):
        sel = SpringShootingSelector(delta_max=1, k_spring=1)
        with pytest.raises(RuntimeError):
            sel.pick(trajectory=self.mytraj)

    def test_impossible_pick(self):
        sel = SpringShootingSelector(delta_max=1, k_spring=1, initial_guess=12)
        sel.pick(trajectory=self.mytraj, direction='forward')
        assert sel._acceptable_snapshot is False

    @property
    def default_selector(self):
        sel = SpringShootingSelector(delta_max=1, k_spring=0,
                                     initial_guess=self.initial_guess)
        sel._fw_prob_list = [1.0, 0.0, 0.0]
        sel._bw_prob_list = [0.0, 0.0, 1.0]
        sel._total_bias = 1.0
        return sel

    def test_f_no_direction(self):
        sel = self.default_selector
        with pytest.raises(NotImplementedError,
                           match="direction"):
            sel.f(self.mytraj[0], self.mytraj)

    def test_f_broken_direction(self):
        sel = self.default_selector
        with pytest.raises(NotImplementedError,
                           match="direction"):
            sel.f(self.mytraj[0], self.mytraj, direction='broken')

    def test_f_no_previous_index(self):
        sel = self.default_selector
        sel.previous_trajectory = self.mytraj
        # Break this to have no initial guess
        sel.previous_snapshot = None
        with pytest.raises(NotImplementedError,
                           match="index"):
            sel.f(self.mytraj[0], self.mytraj, direction='forward')

    def test_f_no_prev_traj(self):
        sel = self.default_selector
        with pytest.raises(NotImplementedError,
                           match="self.previous_trajectory"):
            sel.f(self.mytraj[0], self.mytraj, direction="forward")

    def test_f_wrong_traj(self):
        sel = self.default_selector
        sel.previous_trajectory = self.mytraj
        with pytest.raises(NotImplementedError,
                           match="self.previous_trajectory"):
            # Slice trajectory down to trigger exception
            sel.f(self.mytraj[0], self.mytraj[:4], direction="forward")

    def test_probability_no_direction(self):
        sel = self.default_selector
        with pytest.raises(NotImplementedError,
                           match="direction"):
            sel.probability(self.mytraj[0], self.mytraj)

    def test_probability_broken_direction(self):
        sel = self.default_selector
        with pytest.raises(NotImplementedError,
                           match="direction"):
            sel.probability(self.mytraj[0], self.mytraj, direction='broken')

    def test_probability_forward(self):
        sel = self.default_selector
        sel.previous_trajectory = self.mytraj
        sel._fw_prob_list = [1.0, 1.0, 0.0]
        sel._bw_prob_list = [0.0, 1.0, 1.0]
        sel._fw_total_bias = 2.0
        sel._total_bias = 2.0
        correct = [0.0, 0.0, 0.5, 0.5, 0.0]
        for frame, c in zip(self.mytraj, correct):
            prob = sel.probability(frame, self.mytraj, 'forward')
            assert prob == c

    def test_probability_backward(self):
        sel = self.default_selector
        sel.previous_trajectory = self.mytraj
        sel._fw_prob_list = [1.0, 1.0, 0.0]
        sel._bw_prob_list = [0.0, 1.0, 1.0]
        sel._total_bias = 2.0
        sel._bw_total_bias = 2.0
        # Set index from 3 to min2 as would normally happen
        sel.previous_snapshot = -2
        correct = [0.0, 0.0, 0.0, 0.5, 0.5]
        for frame, c in zip(self.mytraj, correct):
            prob = sel.probability(frame, self.mytraj, 'backward')
            assert prob == c

    def test_forward_pick(self):
        sel = self.default_selector
        pick = sel.pick(trajectory=self.mytraj, direction='forward')
        assert pick == 2
        assert sel.trial_snapshot == 2
        assert sel.previous_snapshot == 3

    def test_backward_pick(self):
        self.initial_guess = 2
        sel = self.default_selector
        pick = sel.pick(trajectory=self.mytraj, direction='backward')
        assert pick == 3
        assert sel.trial_snapshot == -2
        assert sel.previous_snapshot == self.initial_guess

    def test_illegal_forward_pick(self):
        self.initial_guess = 1
        sel = self.default_selector
        pick = sel.pick(trajectory=self.mytraj, direction='forward')
        assert pick == len(self.mytraj)-1
        assert sel.trial_snapshot == len(self.mytraj)-1
        assert sel.previous_snapshot == self.initial_guess
        assert sel._acceptable_snapshot is False

    def test_illegal_backward_pick(self):
        self.initial_guess = 3
        sel = self.default_selector
        pick = sel.pick(trajectory=self.mytraj, direction='backward')
        assert pick == 0
        assert sel.trial_snapshot == -len(self.mytraj)
        assert sel.previous_snapshot == self.initial_guess
        assert sel._acceptable_snapshot is False

    @staticmethod
    def test_failed_loading():
        sel = SpringShootingSelector(delta_max=1, k_spring=0, initial_guess=3)
        details = Details(foo='bar')
        step = FakeStep(details)
        with pytest.raises(RuntimeError):
            sel.restart_from_step(step)

    @staticmethod
    def test_correct_loading():
        sel = SpringShootingSelector(delta_max=1, k_spring=0, initial_guess=3)
        details = Details(initial_trajectory='foo',
                          last_accepted_shooting_index='bar',
                          shooting_index=13,
                          direction='forward')
        step = FakeStep(details)
        sel.restart_from_step(step)
        assert sel.previous_trajectory == 'foo'
        assert sel.previous_snapshot == 'bar'
        assert sel.trial_snapshot == 13
        details.direction = 'backward'
        step = FakeStep(details)
        sel.restart_from_step(step)
        assert sel.trial_snapshot == 10


class MoverTest(SelectorTest):
    def setup(self):
        super(MoverTest, self).setup()
        sel = SpringShootingSelector(delta_max=1, k_spring=0)
        sel._fw_prob_list = [1.0, 0.0, 0.0]
        sel._bw_prob_list = [0.0, 0.0, 1.0]
        sel._total_bias = 1.0
        self.sel = sel
        self.samp = self.ges.samples[0]


class TestSpringMover(MoverTest):

    def test_directionless_call(self):
        mover = SpringMover(ensemble=self.ens, selector=self.sel,
                            engine=self.dyn)

        with pytest.raises(NotImplementedError):
            mover(self.samp)

    def test_stop_before_dynamics(self):
        sel = self.sel
        sel.trial_snapshot = -13
        mover = ForwardSpringMover(ensemble=self.ens, selector=self.sel,
                                   engine=None)
        _, details = mover(self.samp)
        assert details['rejection_reason'] == 'invalid_index'
        assert details['last_accepted_shooting_index'] == -8

    def test_engine_nan_error(self):
        dyn = NaNEngine(None)
        mover = ForwardSpringMover(ensemble=self.ens, selector=self.sel,
                                   engine=dyn)
        with pytest.raises(SampleNaNError):
            _ = mover(self.samp)

    def test_engine_max_length_error(self):
        self.dyn.options['n_frames_max'] = 1
        mover = ForwardSpringMover(ensemble=self.ens, selector=self.sel,
                                   engine=self.dyn)
        with pytest.raises(SampleMaxLengthError):
            _ = mover(self.samp)

    def test_forward_mover(self):
        mover = ForwardSpringMover(ensemble=self.ens, selector=self.sel,
                                   engine=self.dyn)
        trials, details = mover(self.samp)
        traj = trials[0].trajectory
        assert_items_equal([-0.5, 0.1, 0.2, 0.3, 0.4],
                           [s.coordinates[0][0] for s in traj])
        assert details['shooting_index'] == 1
        assert details['last_accepted_shooting_index'] is None
        assert details['direction'] == 'forward'

    def test_backward_mover(self):
        mover = BackwardSpringMover(ensemble=self.ens, selector=self.sel,
                                    engine=self.dyn)
        trials, details = mover(self.samp)
        traj = trials[0].trajectory
        assert_items_equal([-0.1, 0.1, 0.2, 0.3, 0.5],
                           [s.coordinates[0][0] for s in traj])
        assert details['shooting_index'] == 3
        assert details['last_accepted_shooting_index'] is None
        assert details['direction'] == 'backward'


class TestSpringShootingMover(MoverTest):
    def test_correct_init(self):
        mover = SpringShootingMover(ensemble=self.ens, delta_max=1,
                                    k_spring=0.0, engine=self.dyn,
                                    initial_guess=3)
        assert len(mover.movers) == 2

        assert mover.movers[0].selector is mover.movers[1].selector
        assert mover.selector.delta_max == 1
        assert mover.selector.k_spring == 0.0
        assert mover.selector.trial_snapshot == 3
        assert mover.ensemble is self.ens

    def test_from_dict(self):
        dct = {}
        bwmover = BackwardSpringMover(ensemble=self.ens, selector=self.sel,
                                      engine=self.dyn)
        fwmover = ForwardSpringMover(ensemble=self.ens, selector=self.sel,
                                     engine=self.dyn)
        dct['movers'] = [fwmover, bwmover]
        mover = SpringShootingMover.from_dict(dct)
        assert isinstance(mover, SpringShootingMover)
        assert len(mover.movers) == 2
        assert isinstance(mover.movers[0], SpringMover)
        assert isinstance(mover.movers[1], SpringMover)

    def test_dict_cycle(self):
        old_mover = SpringShootingMover(ensemble=self.ens, delta_max=1,
                                        k_spring=0.0, engine=self.dyn,
                                        initial_guess=3)
        dct = old_mover.to_dict()
        mover = SpringShootingMover.from_dict(dct)
        assert isinstance(mover, SpringShootingMover)
        assert len(mover.movers) == 2
        assert isinstance(mover.movers[0], SpringMover)
        assert isinstance(mover.movers[1], SpringMover)
        assert mover.movers[0].selector is mover.movers[1].selector


class TestSpringShootingStrategy(MoveStrategyTestSetup):
    @staticmethod
    def test_init():
        strategy = SpringShootingStrategy(delta_max=1, k_spring=0.0,
                                          initial_guess=3)
        assert strategy.delta_max == 1
        assert strategy.k_spring == 0.0
        assert strategy.initial_guess == 3

    def test_make_movers(self):
        strategy = SpringShootingStrategy(delta_max=1, k_spring=0.0,
                                          initial_guess=3)
        scheme = MoveScheme(self.network)
        shooters = strategy.make_movers(scheme)
        for shooter in shooters:
            assert isinstance(shooter, SpringShootingMover)
            sel = shooter.selector
            assert sel.delta_max == 1
            assert sel.k_spring == 0.0
            assert sel.trial_snapshot == 3


class TestSpringShootingMoveScheme(MoveStrategyTestSetup):
    def test_init(self):
        scheme = SpringShootingMoveScheme(network=self.network,
                                          delta_max=1,
                                          k_spring=0.0,
                                          initial_guess=3)
        assert scheme.delta_max == 1
        assert scheme.k_spring == 0.0
        assert scheme.initial_guess == 3
        for mover in scheme.movers:
            assert isinstance(mover, SpringShootingMover)

    def test_dict_cycle(self):
        scheme = SpringShootingMoveScheme(network=self.network,
                                          delta_max=1,
                                          k_spring=0.0,
                                          initial_guess=3)
        dct = scheme.to_dict()
        scheme2 = SpringShootingMoveScheme.from_dict(dct)
        assert_items_equal(scheme2.movers, scheme.movers)
        assert scheme2.delta_max == 1
        assert scheme2.k_spring == 0.0
        assert scheme2.initial_guess == 3
