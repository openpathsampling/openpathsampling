from builtins import object
from nose.tools import (assert_equal, assert_not_equal, assert_almost_equal,
                        raises)
from nose.plugins.skip import Skip, SkipTest
from openpathsampling.tests.test_helpers import (assert_equal_array_array,
                          assert_not_equal_array_array,
                          make_1d_traj, assert_items_equal,
                          CalvinistDynamics
                         )
from openpathsampling.shooting import *
from openpathsampling.pathmover import ForwardShootMover, BackwardShootMover, SampleMover
from openpathsampling.ensemble import LengthEnsemble
from openpathsampling.sample import Sample, SampleSet

class SelectorTest(object):
    def setup(self):
        self.mytraj = make_1d_traj(coordinates=[-0.5, 0.1, 0.2, 0.3, 0.5],
                                   velocities=[1.0, 1.0, 1.0, 1.0, 1.0])
        self.dyn = CalvinistDynamics([-0.5, -0.4, -0.3, -0.2, -0.1,
                                      0.1, 0.2, 0.3, 0.4, 0.5])
                                      #0.5, 0.4, 0.3, 0.2, 0.1])
        # SampleMover.engine = self.dyn
        self.dyn.initialized = True
        self.ens = LengthEnsemble(5)
        self.gs = SampleSet(Sample(
            replica=0, trajectory=self.mytraj, ensemble=self.ens
        ))


class testFirstFrameSelector(SelectorTest):
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
        assert_items_equal([-0.5, -0.4, -0.3, -0.2, -0.1],
                           [s.coordinates[0][0] for s in samples[0].trajectory]
                          )

class testFinalFrameSelector(SelectorTest):
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
        assert_items_equal([0.1, 0.2, 0.3, 0.4, 0.5],
                           [s.coordinates[0][0] for s in samples[0].trajectory]
                          )
