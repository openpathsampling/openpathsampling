import os
import numpy as np

from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        assert_almost_equal, raises)
from nose.plugins.skip import Skip, SkipTest
from test_helpers import (
    true_func, assert_equal_array_array, make_1d_traj, data_filename
)

import openpathsampling as paths
import logging

logging.getLogger('openpathsampling.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.ensemble').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.netcdfplus').setLevel(logging.CRITICAL)

from openpathsampling.analysis.path_histogram import *
from collections import Counter

class testPathHistogram(object):
    def setup(self):
        self.trajectory = [(0.1, 0.3), (2.1, 3.1), (1.7, 1.4), 
                           (1.6, 0.6), (0.1, 1.4), (2.2, 3.3)]
        self.diag = [(0.25, 0.25), (2.25, 2.25)]

    def test_nointerp_nopertraj(self):
        hist = PathHistogram(left_bin_edges=(0.0, 0.0), 
                             bin_widths=(0.5, 0.5),
                             interpolate=False, per_traj=False)
        hist.add_trajectory(self.trajectory)
        for val in [(0,0), (0,2), (3,1), (3,2)]:
            assert_equal(hist._histogram[val], 1.0)
        assert_equal(hist._histogram[(4,6)], 2.0)
        for val in [(1,0), (1,1), (1,6)]:
            assert_equal(hist._histogram[val], 0.0)

    def test_interp_nopertraj(self):
        hist = PathHistogram(left_bin_edges=(0.0, 0.0), 
                             bin_widths=(0.5, 0.5),
                             interpolate=True, per_traj=False)
        hist.add_trajectory(self.trajectory)
        for val in [(0,0), (0,2), (3,1), (3,2)]:
            assert_equal(hist._histogram[val], 1.0)
        for val in [(0,1), (0,3), (1,4), (2,1), (2,3), (2,5), (3,1), (3,2),
                    (3,3), (3,6)]:
            assert_equal(hist._histogram[val], 1.0)
        for val in [(1,1), (1,2), (1,3), (2,4), (3,4), (4,5), (4,6)]:
            assert_equal(hist._histogram[val], 2.0)
        assert_equal(hist._histogram[(3,5)], 3.0)

    def test_interp_pertraj(self):
        hist = PathHistogram(left_bin_edges=(0.0, 0.0), 
                             bin_widths=(0.5, 0.5),
                             interpolate=True, per_traj=True)
        hist.add_trajectory(self.trajectory)
        for val in [(0,0), (0,2), (3,1), (3,2), (0,1), (0,3), (1,4), (2,1),
                    (2,3), (2,5), (3,1), (3,2), (3,3), (3,6)]:
            assert_equal(hist._histogram[val], 1.0)
        for val in [(1,1), (1,2), (1,3), (2,4), (3,4), (4,5), (4,6)]:
            assert_equal(hist._histogram[val], 1.0)
        assert_equal(hist._histogram[(3,5)], 1.0)

    def test_nointerp_pertraj(self):
        hist = PathHistogram(left_bin_edges=(0.0, 0.0), 
                             bin_widths=(0.5, 0.5),
                             interpolate=False, per_traj=True)
        hist.add_trajectory(self.trajectory)
        for val in [(0,0), (0,2), (3,1), (3,2), (4,6)]:
            assert_equal(hist._histogram[val], 1.0)
        for val in [(1,0), (1,1), (1,6)]:
            assert_equal(hist._histogram[val], 0.0)

    def test_diag_interp(self):
        raise SkipTest

    def test_add_with_weight(self):
        raise SkipTest

