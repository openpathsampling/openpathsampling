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
from opentis.trajectory import Trajectory, Sample

from opentis.ensemble import LengthEnsemble

class testSampleSet(object):
    def setup(self):
        self.ensA = LengthEnsemble(1)
        self.ensB = LengthEnsemble(2)
        traj0A = Trajectory([0.5])
        traj1A = Trajectory([1.0])
        traj2B = Trajectory([0.5, 0.75])
        traj2B_ = Trajectory([0.8, 0.9])
        self.s0A = Sample(replica=0, trajectory=traj0A, ensemble=self.ensA)
        self.s1A = Sample(replica=1, trajectory=traj1A, ensemble=self.ensA)
        self.s2B = Sample(replica=2, trajectory=traj2B, ensemble=self.ensB)
        self.s2B_ = Sample(replica=2, trajectory=traj2B_, ensemble=self.ensB)
        self.testset = SampleSet([self.s0A, self.s1A, self.s2B])

    def test_initialization(self):
        self.testset.consistency_check()

    def test_iter(self):
        samps = [self.s0A, self.s1A, self.s2B]
        for samp in self.testset:
            assert_equal(samp in samps, True)

    def test_len(self):
        assert_equal(len(self.testset), 3)

    def test_getitem_ensemble(self):
        assert_equal(self.testset[self.ensB], self.s2B)
        # TODO: add test that we pick at random for ensA

    def test_getitem_replica(self):
        assert_equal(self.testset[0], self.s0A)
        assert_equal(self.testset[1], self.s1A)
        assert_equal(self.testset[2], self.s2B)
        # TODO: add test that we pick at random

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


