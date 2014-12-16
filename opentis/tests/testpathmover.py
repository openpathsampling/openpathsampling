'''
@author: David W.H. Swenson
'''

import os
import numpy as np

from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        assert_almost_equal, raises)
from nose.plugins.skip import Skip, SkipTest

from opentis.pathmover import *

class testMakeListOfPairs(object):
    def setup(self):
        pass

    def test_not_iterable_type_error(self):
        raise SkipTest
    
    def test_list_of_list_pairs(object):
        raise SkipTest

    def test_list_of_list_notpairs(object):
        raise SkipTest

    def test_list_not_even(object):
        raise SkipTest

    def test_list_even(object):
        raise SkipTest

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

