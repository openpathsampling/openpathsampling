from builtins import object
import collections
from nose.tools import assert_equal, assert_almost_equal, raises
from nose.plugins.skip import Skip, SkipTest
from openpathsampling.tests.test_helpers import (
    assert_equal_array_array, assert_not_equal_array_array, make_1d_traj,
    assert_items_equal, CalvinistDynamics
)
import pytest

from openpathsampling.shooting import *
from openpathsampling.pathmover import ForwardShootMover, BackwardShootMover
from openpathsampling.ensemble import LengthEnsemble
from openpathsampling.sample import Sample, SampleSet
import openpathsampling as paths


class SelectorTest(object):
    def setup(self):
        self.mytraj = make_1d_traj(coordinates=[-0.5, 0.1, 0.2, 0.3, 0.5],
                                   velocities=[1.0, 1.0, 1.0, 1.0, 1.0])
        self.dyn = CalvinistDynamics([-0.5, -0.4, -0.3, -0.2, -0.1,
                                      0.1, 0.2, 0.3, 0.4, 0.5])
        self.dyn.initialized = True
        self.ens = LengthEnsemble(5)
        self.gs = SampleSet(Sample(
            replica=0, trajectory=self.mytraj, ensemble=self.ens
        ))


class TestUniformSelector(SelectorTest):
    def test_pick(self):
        sel = UniformSelector()  # pad_start = pad_end = 1
        pick_idxs = [sel.pick(self.mytraj) for _ in range(100)]
        assert_equal(set(collections.Counter(pick_idxs).keys()), {1, 2, 3})


class TestShootingPointSelector(SelectorTest):
    def test_probability(self):
        # tests this does, in fact, represent the same this as a uniform
        # shooter
        uniform = UniformSelector(pad_start=0, pad_end=0)
        sel = ShootingPointSelector()
        for frame in self.mytraj:
            assert sel.probability(frame, self.mytraj) == \
                    uniform.probability(frame, self.mytraj)

class TestGaussianBiasSelector(SelectorTest):
    def setup(self):
        super(TestGaussianBiasSelector, self).setup()
        self.cv = paths.FunctionCV("Id", lambda x: x.xyz[0][0])
        self.sel = GaussianBiasSelector(self.cv, alpha=2.0, l_0=0.25)
        self.f = [
            0.32465246735834974,  # = exp(-2.0*(-0.5-0.25)**2)
            0.9559974818331,      # = exp(-2.0*(0.1-0.25)**2)
            0.9950124791926823,   # = exp(-2.0*(0.2-0.25)**2)
            0.9950124791926823,   # = exp(-2.0*(0.3-0.25)**2)
            0.8824969025845955,   # = exp(-2.0*(0.5-0.25)**2)
        ]

    def test_pick(self):
        picks = [self.sel.pick(self.mytraj) for _ in range(100)]
        pick_counter = collections.Counter(picks)
        assert set(pick_counter.keys()) == set(range(len(self.mytraj)))
        # final test: 99.5 should happen more than 32.4
        assert pick_counter[2] > pick_counter[0]

    @pytest.mark.parametrize('frame', [0, 1, 2, 3, 4])
    def test_f(self, frame):
        traj = self.mytraj
        assert self.sel.f(traj[frame], traj) == pytest.approx(self.f[frame])

    @pytest.mark.parametrize('frame', [0, 1, 2, 3, 4])
    def test_probability(self, frame):
        norm = sum(self.f)
        traj = self.mytraj
        expected = pytest.approx(self.f[frame] / norm)
        assert self.sel.probability(traj[frame], traj) == expected


class TestFirstFrameSelector(SelectorTest):
    def test_pick(self):
        sel = FirstFrameSelector()
        sp = sel.pick(self.mytraj)
        assert_equal(sp, 0)
        snap = self.mytraj[sp]
        assert_equal(snap.coordinates[0][0], -0.5)

    def test_shooting_move(self):
        self.shooter = ForwardShootMover(
            ensemble=self.ens,
            selector=FirstFrameSelector(),
            engine=self.dyn
        )
        change = self.shooter.move(self.gs)
        samples = change.trials
        assert_equal(len(samples), 1)
        assert_equal(change.accepted, True)
        assert_items_equal(
            [-0.5, -0.4, -0.3, -0.2, -0.1],
            [s.coordinates[0][0] for s in samples[0].trajectory]
        )

    def test_f(self):
        sel = FirstFrameSelector()
        assert sel.f(self.mytraj[0], self.mytraj) == 1
        for frame in self.mytraj[1:]:
            assert sel.f(frame, self.mytraj) == 0


class TestFinalFrameSelector(SelectorTest):
    def test_pick(self):
        sel = FinalFrameSelector()
        sp = sel.pick(self.mytraj)
        assert_equal(sp, 4)
        snap = self.mytraj[sp]
        assert_equal(snap.coordinates[0][0], 0.5)

    def test_shooting_move(self):
        self.shooter = BackwardShootMover(
            ensemble=self.ens,
            selector=FinalFrameSelector(),
            engine=self.dyn
        )
        change = self.shooter.move(self.gs)
        samples = change.trials
        assert_equal(change.accepted, True)
        assert_items_equal(
            [0.1, 0.2, 0.3, 0.4, 0.5],
            [s.coordinates[0][0] for s in samples[0].trajectory]
        )

    def test_f(self):
        sel = FinalFrameSelector()
        for frame in self.mytraj[:-1]:
            assert sel.f(frame, self.mytraj) == 0
        assert sel.f(self.mytraj[-1], self.mytraj) == 1


class TestConstrainedSelector(SelectorTest):
    def setup(self):
        cvx = paths.FunctionCV('ID', lambda snap: snap.xyz[0][0])
        vol = paths.CVDefinedVolume(cvx, float('-inf'), 0)
        self.sel = InterfaceConstrainedSelector(vol)

    def test_pick(self):
        mytraj = make_1d_traj(coordinates=[-0.5, -0.4, -0.3, -0.1,
                                           0.1, 0.2, 0.3, 0.5])
        assert_equal(self.sel.pick(mytraj), 4)

    @raises(RuntimeError)
    def test_all_in_interface(self):
        mytraj = make_1d_traj(coordinates=[-0.5, -0.4, -0.3, -0.1,
                                           -0.1, -0.2, -0.3])
        self.sel.pick(mytraj)

    def test_sum_bias(self):
        mytraj = make_1d_traj(coordinates=[-0.5, 0.1, 0.2, 0.3, 0.5])
        assert_equal(self.sel.sum_bias(mytraj), 1.0)

    def test_f(self):
        mytraj = make_1d_traj(coordinates=[-0.5, -0.4, -0.3, -0.1,
                                           0.1, 0.2, 0.3, 0.5])
        expected_idx = 4
        idx = self.sel.pick(mytraj)
        frame = mytraj[idx]
        assert_equal(self.sel.f(frame, mytraj), 1.0)
        for idx1, frame in enumerate(mytraj):
            if (idx1 != expected_idx):
                assert_equal(self.sel.f(frame, mytraj), 0.0)
