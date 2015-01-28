
from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        assert_almost_equal, raises)
from nose.plugins.skip import Skip, SkipTest
from test_helpers import (assert_equal_array_array, 
                          assert_not_equal_array_array,
                          make_1d_traj
                         )
from opentis.shooting import *


class testFirstFrameSelector(object):
    def setup(self):
        pass
    
    def test_pick(self):
        raise SkipTest

    def test_shooting_move(self):
        raise SkipTest

class testFinalFrameSelector(object):
    def setup(self):
        pass

    def test_pick(self):
        raise SkipTest

    def test_shooting_move(self):
        raise SkipTest
