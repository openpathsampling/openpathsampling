
from builtins import object
from nose.tools import assert_equal, assert_not_equal, raises
from nose.plugins.skip import Skip, SkipTest
from openpathsampling.range_logic import *

class TestRangeLogic(object):
    def test_range_and(self):
        assert_equal(range_and(1, 3, 2, 4), [(2, 3)])
        assert_equal(range_and(2, 4, 1, 3), [(2, 3)])
        assert_equal(range_and(1, 2, 3, 4), None)
        assert_equal(range_and(3, 4, 1, 2), None)
        assert_equal(range_and(1, 4, 2, 3), [(2, 3)])
        assert_equal(range_and(2, 3, 1, 4), [(2, 3)])
        assert_equal(range_and(1, 2, 1, 2), 1)

    def test_range_or(self):
        assert_equal(range_or(1, 3, 2, 4), [(1, 4)])
        assert_equal(range_or(2, 4, 1, 3), [(1, 4)])
        assert_equal(range_or(1, 2, 3, 4), [(1, 2), (3, 4)])
        assert_equal(range_or(3, 4, 1, 2), [(3, 4), (1, 2)])
        assert_equal(range_or(1, 4, 2, 3), [(1, 4)])
        assert_equal(range_or(2, 3, 1, 4), [(1, 4)])
        assert_equal(range_or(1, 2, 1, 2), 1)


    def test_range_sub(self):
        assert_equal(range_sub(1, 3, 2, 4), [(1, 2)])
        assert_equal(range_sub(2, 4, 1, 3), [(3, 4)])
        assert_equal(range_sub(1, 2, 3, 4), 1)
        assert_equal(range_sub(3, 4, 1, 2), 1)
        assert_equal(range_sub(1, 4, 2, 3), [(1, 2), (3, 4)])
        assert_equal(range_sub(2, 3, 1, 4), None)
        assert_equal(range_sub(1, 2, 1, 2), None)
        assert_equal(range_sub(0.1, 0.4, 0.1, 0.3), [(0.3, 0.4)])


class TestPeriodicRangeLogic(object):

    def test_periodic_order(self):
        # orders without wrapping
        assert_equal(periodic_ordering(1, 2, 3, 4), [0, 1, 2, 3])
        assert_equal(periodic_ordering(1, 3, 2, 4), [0, 2, 1, 3])
        assert_equal(periodic_ordering(4, 3, 2, 1), [0, 3, 2, 1])
        assert_equal(periodic_ordering(1, 2, 1, 2), [0, 2, 1, 3])
        assert_equal(periodic_ordering(2, 4, 1, 3), [1, 3, 0, 2])
        assert_equal(periodic_ordering(1, 2, 4, 3), [1, 2, 0, 3])

    def test_periodic_and(self):
        assert_equal(periodic_range_and(0.1, 0.3, 0.2, 0.4), [(0.2, 0.3)])
        assert_equal(periodic_range_and(0.2, 0.4, 0.1, 0.3), [(0.2, 0.3)])
        assert_equal(periodic_range_and(1, 2, 3, 4), None)
        assert_equal(periodic_range_and(3, 4, 1, 2), None)
        assert_equal(periodic_range_and(1, 4, 2, 3), [(2, 3)])
        assert_equal(periodic_range_and(2, 3, 1, 4), [(2, 3)])
        assert_equal(periodic_range_and(1, 2, 1, 2), 1)
        assert_equal(periodic_range_and(1, 2, 2, 1), None)
        assert_equal(periodic_range_and(2, 1, 1, 4), [(2, 4)])
        assert_equal(periodic_range_and(0.1, 0.4, 0.3, 0.2), 
                     [(0.1, 0.2), (0.3, 0.4)])

    def test_periodic_or(self):
        assert_equal(periodic_range_or(0.1, 0.3, 0.2, 0.4), [(0.1, 0.4)])
        assert_equal(periodic_range_or(0.2, 0.4, 0.1, 0.3), [(0.1, 0.4)])
        assert_equal(periodic_range_or(1, 2, 3, 4), [(1, 2), (3, 4)])
        assert_equal(periodic_range_or(3, 4, 1, 2), [(3, 4), (1, 2)])
        assert_equal(periodic_range_or(1, 4, 2, 3), [(1, 4)])
        assert_equal(periodic_range_or(2, 3, 1, 4), [(1, 4)])
        assert_equal(periodic_range_or(1, 2, 1, 2), 1)
        assert_equal(periodic_range_or(1, 2, 2, 1), -1)
        assert_equal(periodic_range_or(0.1, 0.4, 0.3, 0.2), -1)
        assert_equal(periodic_range_or(2, 1, 1, 4), -1)


    def test_periodic_sub(self):
        assert_equal(periodic_range_sub(0.1, 0.3, 0.2, 0.4), [(0.1, 0.2)])
        assert_equal(periodic_range_sub(0.2, 0.4, 0.1, 0.3), [(0.3, 0.4)])
        assert_equal(periodic_range_sub(1, 2, 3, 4), 1)
        assert_equal(periodic_range_sub(3, 4, 1, 2), 1)
        assert_equal(periodic_range_sub(1, 4, 2, 3), [(1, 2), (3, 4)])
        assert_equal(periodic_range_sub(2, 3, 1, 4), None)
        assert_equal(periodic_range_sub(1, 2, 1, 2), None)
        assert_equal(periodic_range_sub(1, 2, 2, 1), 1)
        assert_equal(periodic_range_sub(2, 1, 1, 4), [(4, 1)])
        assert_equal(periodic_range_sub(0.1, 0.4, 0.3, 0.2), [(0.2, 0.3)])
        assert_equal(periodic_range_sub(0.1, 0.4, 0.1, 0.3), [(0.3, 0.4)])
