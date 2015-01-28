
from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        assert_almost_equal, raises)
from nose.plugins.skip import Skip, SkipTest
from test_helpers import (assert_equal_array_array, 
                          assert_not_equal_array_array,
                          make_1d_traj
                         )
from opentis.shooting import *

from opentis.pathmover import PathMover, OneWayShootingMover

class SelectorTest(object):
    def setup(self):
        self.mytraj = make_1d_traj([-0.5, 0.1, 0.3, 0.5])
        #self.shooter = OneWayShootingMover(FirstFrameSelector())


class testFirstFrameSelector(SelectorTest):
    def test_pick(self):
        sel = FirstFrameSelector()
        sp = sel.pick(self.mytraj)
        assert_equal(sp.selector, sel)
        assert_equal(sp.trajectory, self.mytraj)
        assert_equal(sp.index, 0)
        assert_equal(sp.f, 1.0)
        assert_equal(sp.sum_bias, 1.0)
        snap = sp.snapshot
        assert_equal(snap.coordinates[0][0], -0.5)

    def test_shooting_move(self):
        # this just tests whether the shooting move runs, not correctness
        raise SkipTest

class testFinalFrameSelector(SelectorTest):
    def test_pick(self):
        raise SkipTest

    def test_shooting_move(self):
        raise SkipTest
