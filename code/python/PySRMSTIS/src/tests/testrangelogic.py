
import os
import sys
from nose.tools import assert_equal, assert_not_equal, raises
from nose.plugins.skip import Skip, SkipTest
sys.path.append(os.path.abspath('../'))
from range_logic import *

class testRangeLogic(object):
    def test_range_and(self):
        assert_equal(range_and(1, 3, 2, 4), (2, 3))
        assert_equal(range_and(2, 4, 1, 3), (2, 3))
        assert_equal(range_and(1, 2, 3, 4), None)
        assert_equal(range_and(3, 4, 1, 2), None)
        assert_equal(range_and(1, 4, 2, 3), (2, 3))
        assert_equal(range_and(2, 3, 1, 4), (2, 3))
        assert_equal(range_and(1, 2, 1, 2), 1)

    def test_range_or(self):
        assert_equal(range_or(1, 3, 2, 4), (1, 4))
        assert_equal(range_or(2, 4, 1, 3), (1, 4))
        assert_equal(range_or(1, 2, 3, 4), 2)
        assert_equal(range_or(3, 4, 1, 2), 2)
        assert_equal(range_or(1, 4, 2, 3), (1, 4))
        assert_equal(range_or(2, 3, 1, 4), (1, 4))
        assert_equal(range_or(1, 2, 1, 2), 1)


    def test_range_sub(self):
        assert_equal(range_sub(1, 3, 2, 4), (1, 2))
        assert_equal(range_sub(2, 4, 1, 3), (3, 4))
        assert_equal(range_sub(1, 2, 3, 4), 1)
        assert_equal(range_sub(3, 4, 1, 2), 1)
        assert_equal(range_sub(1, 4, 2, 3), 2)
        assert_equal(range_sub(2, 3, 1, 4), None)
        assert_equal(range_sub(1, 2, 1, 2), None)

