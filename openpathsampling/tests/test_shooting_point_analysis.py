from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        raises, assert_almost_equal, assert_true)
from nose.plugins.skip import SkipTest
from numpy.testing import assert_array_almost_equal
from test_helpers import make_1d_traj

import openpathsampling as paths
import openpathsampling.engines as peng
import numpy as np

from openpathsampling.analysis.shooting_point_analysis import *

class testTransformedDict(object):
    def setup(self):
        self.untransformed = {(0, 1) : "a", (1, 2) : "b", (2, 3) : "c"}
        self.transformed = {0 : "a", 1 : "b", 2 : "c"}
        self.hash_function = lambda x : x[0]

    def test_initialization(self):
        test_dict = TransformedDict(self.hash_function, self.untransformed)
        assert_equal(test_dict.store, self.transformed)
        assert_equal(test_dict.hash_representatives, 
                     {0: (0,1), 1: (1,2), 2: (2,3)})

    def test_set_get(self):
        pass

    def test_update(self):
        pass

    def test_del(self):
        pass

    def test_iter(self):
        pass

    def test_len(self):
        pass

