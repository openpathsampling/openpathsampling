from __future__ import absolute_import
from builtins import object
import os
import numpy as np
import collections

import pytest

from nose.tools import (assert_equal, assert_not_equal, assert_almost_equal,
                        raises)
from nose.plugins.skip import Skip, SkipTest
from .test_helpers import (
    true_func, assert_equal_array_array, make_1d_traj, data_filename,
    assert_items_almost_equal
)

import openpathsampling as paths
import logging

logging.getLogger('openpathsampling.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.ensemble').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.netcdfplus').setLevel(logging.CRITICAL)

from openpathsampling.analysis.path_histogram import *
from collections import Counter

class PathHistogramTester(object):
    def setup(self):
        self.HAS_TQDM = paths.progress.HAS_TQDM
        paths.progress.HAS_TQDM = False
        self.trajectory = [(0.1, 0.3), (2.1, 3.1), (1.7, 1.4),
                           (1.6, 0.6), (0.1, 1.4), (2.2, 3.3)]
        self.diag = [(0.25, 0.25), (2.25, 2.25)]

    def teardown(self):
        paths.progress.HAS_TQDM = self.HAS_TQDM

    def test_nopertraj(self):
        counter = collections.Counter(self.expected_bins)
        hist = PathHistogram(left_bin_edges=(0.0, 0.0),
                             bin_widths=(0.5, 0.5),
                             interpolate=self.Interpolator, per_traj=False)
        hist.add_trajectory(self.trajectory)
        assert hist._histogram == counter

    def test_pertraj(self):
        counter = collections.Counter(set(self.expected_bins))
        hist = PathHistogram(left_bin_edges=(0.0, 0.0),
                             bin_widths=(0.5, 0.5),
                             interpolate=self.Interpolator, per_traj=True)
        hist.add_trajectory(self.trajectory)
        assert hist._histogram == counter

    def test_diag(self):
        hist = PathHistogram(left_bin_edges=(0.0, 0.0),
                             bin_widths=(0.5, 0.5),
                             interpolate=self.Interpolator, per_traj=True)
        hist.add_trajectory(self.diag)
        expected = collections.Counter([(i, i) for i in range(5)])
        assert hist._histogram == expected

    def test_same_cell(self):
        # check interpolation if successive frames in same cell
        traj = [(0.1, 0.3), (0.2, 0.2), (0.4, 0.6), (0.3, 0.1)]
        hist = PathHistogram(left_bin_edges=(0.0, 0.0),
                             bin_widths=(0.5, 0.5),
                             interpolate=self.Interpolator, per_traj=False)
        hist.add_trajectory(traj)
        assert_equal(len(list(hist._histogram.keys())), 2)
        assert_equal(hist._histogram[(0,0)], 3)
        assert_equal(hist._histogram[(0,1)], 1)


class TestPathHistogramNoInterpolate(PathHistogramTester):
    Interpolator = NoInterpolation
    def setup(self):
        super(TestPathHistogramNoInterpolate, self).setup()
        self.expected_bins = [(0, 0), (4, 6), (3, 2),
                              (3, 1), (0, 2), (4, 6)]

    def test_diag(self):
        hist = PathHistogram(left_bin_edges=(0.0, 0.0),
                             bin_widths=(0.5, 0.5),
                             interpolate=self.Interpolator, per_traj=True)
        hist.add_trajectory(self.diag)
        assert hist._histogram == collections.Counter([(0, 0), (4, 4)])


class TestPathHistogramSubdivideInterpolate(PathHistogramTester):
    Interpolator = SubdivideInterpolation
    def setup(self):
        super(TestPathHistogramSubdivideInterpolate, self).setup()
        self.expected_bins = [
            (0, 0),  # initial
            (0, 1), (1, 1), (1, 2), (1, 3), (2, 3), (2, 4), (3, 4),
            (3, 5), (4, 5), (4, 6),  # 0->1
            (4, 5), (3, 5), (3, 4), (3, 3), (3, 2),  # 1->2
            (3, 1),  # 2->3
            (2, 1), (1, 1), (1, 2), (0, 2),  # 3->4
            (0, 3), (1, 3), (1, 4), (2, 4), (2, 5), (3, 5),
            (3, 6), (4, 6)  # 4->5
        ]


class TestPathHistogramBesenhamLikeInterpolate(PathHistogramTester):
    Interpolator = BresenhamLikeInterpolation
    def setup(self):
        super(TestPathHistogramBesenhamLikeInterpolate, self).setup()
        # NOTE: the (4, 5) in the 1->2 transition is exactly on bin edge;
        # left bin edge is inclusive
        self.expected_bins = [
            (0, 0),  # initial
            (0, 1), (1, 2), (2, 3), (2, 4), (3, 5), (4, 6),  # 0->1
            (4, 5), (3, 4), (3, 3), (3, 2),  # 1->2
            (3, 1),  # 2->3
            (2, 1), (1, 2), (0, 2),  # 3->4
            (1, 3), (2, 4), (3, 5), (4, 6)  # 4->5
        ]


class TestPathHistogramBesenhamInterpolate(PathHistogramTester):
    Interpolator = BresenhamInterpolation
    def setup(self):
        super(TestPathHistogramBesenhamInterpolate, self).setup()
        self.expected_bins = [
            (0, 0),  # initial
            (1, 1), (1, 2), (2, 3), (3, 4), (3, 5), (4, 6),  # 0->1
            (4, 5), (4, 4), (3, 3), (3, 2),  # 1->2
            (3, 1),  # 2->3
            (2, 1), (1, 2), (0, 2),  # 3->4
            (1, 3), (2, 4), (3, 5), (4, 6)  # 4->5
        ]

class TestPathHistogram(object):
    # tests of fundamental things in PathHistogram, not interpolators
    def setup(self):
        self.HAS_TQDM = paths.progress.HAS_TQDM
        paths.progress.HAS_TQDM = False
        self.trajectory = [(0.1, 0.3), (2.1, 3.1), (1.7, 1.4),
                           (1.6, 0.6), (0.1, 1.4), (2.2, 3.3)]
        self.diag = [(0.25, 0.25), (2.25, 2.25)]

    def teardown(self):
        paths.progress.HAS_TQDM = self.HAS_TQDM

    def test_add_with_weight(self):
        hist = PathHistogram(left_bin_edges=(0.0, 0.0),
                             bin_widths=(0.5, 0.5),
                             interpolate=SubdivideInterpolation,
                             per_traj=True)
        hist.add_trajectory(self.trajectory)
        hist.add_trajectory(self.diag, weight=2)
        for val in [(0,1), (0,2), (0,3), (1,2), (1,3), (1,4), (2,1), (2,3),
                    (2,4), (2,5), (3,1), (3,2), (3,4), (3,5), (3,6), (4,5),
                    (4,6)]:
            assert_equal(hist._histogram[val], 1.0)
        for val in [(2,2), (4,4)]:
            assert_equal(hist._histogram[val], 2.0)
        for val in [(0,0), (1,1), (3,3)]:
            assert_equal(hist._histogram[val], 3.0)
        for val in [(0,4), (0,5), (0.6), (0,7), (-1,0)]:
            assert_equal(hist._histogram[val], 0.0)

    def test_add_data_to_histograms(self):
        hist = PathHistogram(left_bin_edges=(0.0, 0.0),
                             bin_widths=(0.5, 0.5),
                             interpolate=SubdivideInterpolation,
                             per_traj=True)
        counter = hist.add_data_to_histogram([self.trajectory, self.diag],
                                             weights=[1.0, 2.0])
        for val in [(0,1), (0,2), (0,3), (1,2), (1,3), (1,4), (2,1), (2,3),
                    (2,4), (2,5), (3,1), (3,2), (3,4), (3,5), (3,6), (4,5),
                    (4,6)]:
            assert_equal(counter[val], 1.0)
        for val in [(2,2), (4,4)]:
            assert_equal(counter[val], 2.0)
        for val in [(0,0), (1,1), (3,3)]:
            assert_equal(counter[val], 3.0)
        for val in [(0,4), (0,5), (0.6), (0,7), (-1,0)]:
            assert_equal(counter[val], 0.0)

    def test_add_data_to_histograms_no_weight(self):
        hist = PathHistogram(left_bin_edges=(0.0, 0.0),
                             bin_widths=(0.5, 0.5),
                             interpolate=SubdivideInterpolation,
                             per_traj=True)
        counter = hist.add_data_to_histogram([self.trajectory, self.diag])
        for val in [(0,1), (0,2), (0,3), (1,2), (1,3), (1,4), (2,1), (2,3),
                    (2,4), (2,5), (3,1), (3,2), (3,4), (3,5), (3,6), (4,5),
                    (4,6)]:
            assert_equal(counter[val], 1.0)
        for val in [(2,2), (4,4)]:
            assert_equal(counter[val], 1.0)
        for val in [(0,0), (1,1), (3,3)]:
            assert_equal(counter[val], 2.0)
        for val in [(0,4), (0,5), (0.6), (0,7), (-1,0)]:
            assert_equal(counter[val], 0.0)


class TestPathDensityHistogram(object):
    def setup(self):
        self.HAS_TQDM = paths.progress.HAS_TQDM
        paths.progress.HAS_TQDM = False
        id_cv = paths.FunctionCV("Id",
                                 lambda snap : snap.xyz[0][0])
        sin_cv = paths.FunctionCV("sin",
                                  lambda snap : np.sin(snap.xyz[0][0]))
        square_cv = paths.FunctionCV("x^2",
                                     lambda snap : snap.xyz[0][0]**2)
        self.cvs = [id_cv, sin_cv, square_cv]
        self.left_bin_edges = (0,0,0)
        self.bin_widths = (0.25, 0.5, 0.4)

        self.traj1 = make_1d_traj([0.1, 0.51, 0.61])
        # [(0, 0, 0), (2, 0, 0), (2, 1, 0)]
        # interpolate: (1, 0, 0)
        self.traj2 = make_1d_traj([0.6, 0.7])
        # [(2, 1, 0), (2, 1, 1)]

    def teardown(self):
        paths.progress.HAS_TQDM = self.HAS_TQDM

    def test_histogram_no_weights(self):
        hist = PathDensityHistogram(self.cvs, self.left_bin_edges,
                                    self.bin_widths)
        counter = hist.histogram([self.traj1, self.traj2])
        assert_equal(len(counter), 5)
        for bin_label in [(0,0,0), (2,0,0), (1,0,0), (2,1,1)]:
            assert_equal(counter[bin_label], 1.0)
        assert_equal(counter[(2,1,0)], 2.0)
        for bin_label in [(1,1,1), (2,2,2), (-1,0,0), (0,-1,0)]:
            assert_equal(counter[bin_label], 0.0)

    def test_histogram_with_weights(self):
        hist = PathDensityHistogram(self.cvs, self.left_bin_edges,
                                    self.bin_widths)
        counter = hist.histogram([self.traj1, self.traj2], 
                                 weights=[1.0, 2.0])
        assert_equal(len(counter), 5)
        for bin_label in [(0,0,0), (2,0,0), (1,0,0)]:
            assert_equal(counter[bin_label], 1.0)
        assert_equal(counter[(2,1,1)], 2.0)
        assert_equal(counter[(2,1,0)], 3.0)
        for bin_label in [(1,1,1), (2,2,2), (-1,0,0), (0,-1,0)]:
            assert_equal(counter[bin_label], 0.0)

    def test_histogram_single_traj(self):
        hist = PathDensityHistogram(self.cvs, self.left_bin_edges,
                                    self.bin_widths)
        counter = hist.histogram(self.traj1)
        assert_equal(len(counter), 4)
        for bin_label in [(0,0,0), (2,0,0), (1,0,0), (2,1,0)]:
            assert_equal(counter[bin_label], 1.0)
        for bin_label in [(1,1,1), (2,2,2), (-1,0,0), (0,-1,0)]:
            assert_equal(counter[bin_label], 0.0)

    def test_map_to_float_bins_trajectory(self):
        hist = PathDensityHistogram(self.cvs, self.left_bin_edges,
                                    self.bin_widths)
        bin_traj = hist.map_to_float_bins(self.traj2)
        # cv_traj = [[0.6, 0.564642473395035, 0.36],
        #            [0.7, 0.644217687237691, 0.49]]
        assert_items_almost_equal(bin_traj[0], [2.4, 1.12928494679007, 0.9])
        assert_items_almost_equal(bin_traj[1], [2.8, 1.28843537447538, 1.225])

    def test_map_to_float_bins_cvtraj(self):
        hist = PathDensityHistogram(self.cvs, self.left_bin_edges,
                                    self.bin_widths)
        cv_traj = [[0.6, 0.564642473395035, 0.36],
                   [0.7, 0.644217687237691, 0.49]]
        bin_traj = hist.map_to_float_bins(cv_traj)
        assert_items_almost_equal(bin_traj[0], [2.4, 1.12928494679007, 0.9])
        assert_items_almost_equal(bin_traj[1], [2.8, 1.28843537447538, 1.225])

