from nose.tools import (assert_equal, assert_not_equal, assert_items_equal, 
                        assert_almost_equal, raises)
from nose.plugins.skip import SkipTest
from test_helpers import assert_items_almost_equal

import logging

from openpathsampling.analysis import LookupFunction, LookupFunctionGroup

class testLookupFunctionGroup(object):
    def setup(self):
        x1 = [0.0, 1.0, 2.0, 3.0, 4.0]
        y1 = [0.1, 1.3, 1.8, 3.3, 4.2]
        luf1 = LookupFunction(x1, y1)
        x2 = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
        y2 = [0.6, 0.8, 1.3, 1.9, 2.6, 2.8, 3.3, 3.9]
        luf2 = LookupFunction(x2, y2)
        x3 = [1.0, 2.0]
        y3 = [1.0, 2.0]
        luf3 = LookupFunction(x3, y3)
        y4 = [1.1, 1.9]
        luf4 = LookupFunction(x3, y4)
        y5 = [1.2, 2.2]
        luf5 = LookupFunction(x3, y5)
        self.group = LookupFunctionGroup([luf1, luf2, luf3, luf4, luf5])

    def test_all_x(self):
        self.group.use_x = "all"
        assert_items_equal(
            self.group.x, [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
        )

    def test_shared_x(self):
        self.group.use_x = "shared"
        assert_items_equal(self.group.x, [1.0, 2.0])

    def test_mean_shared(self):
        self.group.use_x = "shared"
        assert_almost_equal(self.group.mean(1.0), 1.08)
        assert_almost_equal(self.group(1.0), 1.08)

    def test_mean_allx(self):
        self.group.use_x = "all"
        raise SkipTest

    def test_std_shared(self):
        self.group.use_x = "shared"
        assert_almost_equal(self.group.std(1.0), 0.17204650534085253)

    def test_std_allx(self):
        self.group.use_x = "all"
        raise SkipTest

    def test_getitem(self):
        raise SkipTest

    def test_setitem(self):
        raise SkipTest
