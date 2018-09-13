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
import openpathsampling as paths

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
        assert_items_equal([-0.5, -0.4, -0.3, -0.2, -0.1],
                           [s.coordinates[0][0] for s in samples[0].trajectory]
                          )

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
        assert_items_equal([0.1, 0.2, 0.3, 0.4, 0.5],
                           [s.coordinates[0][0] for s in samples[0].trajectory]
                          )

class TestConstrainedSelector(SelectorTest): 
    def setup(self):
        cvx= paths.FunctionCV('ID',lambda snap:snap.xyz[0][0])
        vol= paths.CVDefinedVolume(cvx,float('-inf'),0 )
        self.sel = InterfaceConstrainedSelector(vol)
        
    def test_pick(self):
        mytraj = make_1d_traj(coordinates=[-0.5,-0.4,-0.3,-0.1, 0.1, 0.2, 0.3, 0.5])
        assert_equal(self.sel.pick(mytraj),4)

    @raises(RuntimeError)    
    def test_allininterface(self):
        mytraj = make_1d_traj(coordinates=[-0.5,-0.4,-0.3,-0.1, -0.1, -0.2, -0.3, -0.5])
        self.sel.pick(mytraj)
        
 
    def test_sum_bias(self):
        mytraj = make_1d_traj(coordinates=[-0.5, 0.1, 0.2, 0.3, 0.5])
        assert_equal(self.sel.sum_bias(mytraj),1.0)
        
    def test_f(self):
        mytraj = make_1d_traj(coordinates=[-0.5,-0.4,-0.3,-0.1, 0.1, 0.2, 0.3, 0.5])
        expected_idx=4;
        idx=self.sel.pick(mytraj)
        frame = mytraj[idx]
        assert_equal(self.sel.f(frame,mytraj),1.0)
        for idx1,frame in enumerate(mytraj):
            if (idx1 != expected_idx): 
                assert_equal(self.sel.f(frame,mytraj),0.0)
