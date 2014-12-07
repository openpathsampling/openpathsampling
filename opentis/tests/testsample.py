'''
@author: David W.H. Swenson
'''

import os
import numpy as np

from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        assert_almost_equal, raises)
from nose.plugins.skip import Skip, SkipTest
from test_helpers import true_func, assert_equal_array_array

from opentis.sample import *
from opentis.trajectory import Sample

class testSampleSet(object):
    def setup(self):
        pass

    def test_initialization(self):
        raise SkipTest

    def test_iter(self):
        raise SkipTest

    def test_len(self):
        raise SkipTest

    def test_getitem_ensemble(self):
        raise SkipTest

    def test_getitem_replica(self):
        raise SkipTest

    def test_setitem_ensemble(self):
        raise SkipTest

    def test_setitem_replica(self):
        raise SkipTest

    def test_setitem_itemexists(self):
        raise SkipTest

    def test_additem(self):
        raise SkipTest

    def test_additem_itemexists(self):
        raise SkipTest

    def test_deleteitem(self):
        raise SkipTest

    def test_apply_samples(self):
        raise SkipTest


