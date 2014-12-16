'''
@author: David W.H. Swenson
'''

import os
import numpy as np

from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        assert_almost_equal, raises)
from nose.plugins.skip import Skip, SkipTest
from test_helpers import assert_equal_array_array, assert_not_equal_array_array

from opentis.pathmover import *

class testMakeListOfPairs(object):
    def setup(self):
        self.correct = [ [0, 1], [2, 3], [4, 5] ]
        pass

    @raises(TypeError)
    def test_not_iterable_type_error(self):
        result = make_list_of_pairs(1)
    
    def test_list_of_list_pairs(self):
        result = make_list_of_pairs([[0, 1], [2, 3], [4, 5]])
        assert_equal_array_array(result, self.correct)

    @raises(AssertionError)
    def test_list_of_list_notpairs(self):
        result = make_list_of_pairs([[0], [1], [2]])

    @raises(AssertionError)
    def test_list_not_even(self):
        result = make_list_of_pairs([0, 1, 2, 3, 4])

    def test_list_even(self):
        result = make_list_of_pairs([0, 1, 2, 3, 4, 5])
        assert_equal_array_array(result, self.correct)

class testPathMover(object):
    def test_legal_sample_set(self):
        raise SkipTest

    def test_select_sample(self):
        raise SkipTest

class testForwardShootMover(object):
    def setup(self):
        pass

    def test_move(self):
        raise SkipTest

class testBackwardShootMover(object):
    def setup(self):
        pass

    def test_move(self):
        raise SkipTest

