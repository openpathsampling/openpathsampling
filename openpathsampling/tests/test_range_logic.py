
from builtins import object
from openpathsampling.range_logic import *

class TestRangeLogic(object):
    def test_range_and(self):
        assert range_and(1, 3, 2, 4) == [(2, 3)]
        assert range_and(2, 4, 1, 3) == [(2, 3)]
        assert range_and(1, 2, 3, 4) is None
        assert range_and(3, 4, 1, 2) is None
        assert range_and(1, 4, 2, 3) == [(2, 3)]
        assert range_and(2, 3, 1, 4) == [(2, 3)]
        assert range_and(1, 2, 1, 2) == 1

    def test_range_or(self):
        assert range_or(1, 3, 2, 4) == [(1, 4)]
        assert range_or(2, 4, 1, 3) == [(1, 4)]
        assert range_or(1, 2, 3, 4) == [(1, 2), (3, 4)]
        assert range_or(3, 4, 1, 2) == [(3, 4), (1, 2)]
        assert range_or(1, 4, 2, 3) == [(1, 4)]
        assert range_or(2, 3, 1, 4) == [(1, 4)]
        assert range_or(1, 2, 1, 2) == 1


    def test_range_sub(self):
        assert range_sub(1, 3, 2, 4) == [(1, 2)]
        assert range_sub(2, 4, 1, 3) == [(3, 4)]
        assert range_sub(1, 2, 3, 4) == 1
        assert range_sub(3, 4, 1, 2) == 1
        assert range_sub(1, 4, 2, 3) == [(1, 2), (3, 4)]
        assert range_sub(2, 3, 1, 4) is None
        assert range_sub(1, 2, 1, 2) is None
        assert range_sub(0.1, 0.4, 0.1, 0.3) == [(0.3, 0.4)]


class TestPeriodicRangeLogic(object):

    def test_periodic_order(self):
        # orders without wrapping
        assert periodic_ordering(1, 2, 3, 4) == [0, 1, 2, 3]
        assert periodic_ordering(1, 3, 2, 4) == [0, 2, 1, 3]
        assert periodic_ordering(4, 3, 2, 1) == [0, 3, 2, 1]
        assert periodic_ordering(1, 2, 1, 2) == [0, 2, 1, 3]
        assert periodic_ordering(2, 4, 1, 3) == [1, 3, 0, 2]
        assert periodic_ordering(1, 2, 4, 3) == [1, 2, 0, 3]

    def test_periodic_and(self):
        assert periodic_range_and(0.1, 0.3, 0.2, 0.4) == [(0.2, 0.3)]
        assert periodic_range_and(0.2, 0.4, 0.1, 0.3) == [(0.2, 0.3)]
        assert periodic_range_and(1, 2, 3, 4) is None
        assert periodic_range_and(3, 4, 1, 2) is None
        assert periodic_range_and(1, 4, 2, 3) == [(2, 3)]
        assert periodic_range_and(2, 3, 1, 4) == [(2, 3)]
        assert periodic_range_and(1, 2, 1, 2) == 1
        assert periodic_range_and(1, 2, 2, 1) is None
        assert periodic_range_and(2, 1, 1, 4) == [(2, 4)]
        assert periodic_range_and(0.1, 0.4, 0.3, 0.2) \
                == [(0.1, 0.2), (0.3, 0.4)]

    def test_periodic_or(self):
        assert periodic_range_or(0.1, 0.3, 0.2, 0.4) == [(0.1, 0.4)]
        assert periodic_range_or(0.2, 0.4, 0.1, 0.3) == [(0.1, 0.4)]
        assert periodic_range_or(1, 2, 3, 4) == [(1, 2), (3, 4)]
        assert periodic_range_or(3, 4, 1, 2) == [(3, 4), (1, 2)]
        assert periodic_range_or(1, 4, 2, 3) == [(1, 4)]
        assert periodic_range_or(2, 3, 1, 4) == [(1, 4)]
        assert periodic_range_or(1, 2, 1, 2) == 1
        assert periodic_range_or(1, 2, 2, 1) == -1
        assert periodic_range_or(0.1, 0.4, 0.3, 0.2) == -1
        assert periodic_range_or(2, 1, 1, 4) == -1


    def test_periodic_sub(self):
        assert periodic_range_sub(0.1, 0.3, 0.2, 0.4) == [(0.1, 0.2)]
        assert periodic_range_sub(0.2, 0.4, 0.1, 0.3) == [(0.3, 0.4)]
        assert periodic_range_sub(1, 2, 3, 4) == 1
        assert periodic_range_sub(3, 4, 1, 2) == 1
        assert periodic_range_sub(1, 4, 2, 3) == [(1, 2), (3, 4)]
        assert periodic_range_sub(2, 3, 1, 4) is None
        assert periodic_range_sub(1, 2, 1, 2) is None
        assert periodic_range_sub(1, 2, 2, 1) == 1
        assert periodic_range_sub(2, 1, 1, 4) == [(4, 1)]
        assert periodic_range_sub(0.1, 0.4, 0.3, 0.2) == [(0.2, 0.3)]
        assert periodic_range_sub(0.1, 0.4, 0.1, 0.3) == [(0.3, 0.4)]
