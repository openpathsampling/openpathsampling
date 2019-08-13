from __future__ import division
from __future__ import absolute_import
from past.utils import old_div
from builtins import object
from nose.tools import (assert_equal, assert_not_equal, raises,
                        assert_almost_equal)
from nose.plugins.skip import SkipTest
from .test_helpers import assert_items_almost_equal, assert_items_equal
import pytest
import logging
logging.getLogger('openpathsampling.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.ensemble').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.netcdfplus').setLevel(logging.CRITICAL)

import collections

from openpathsampling.numerics import (Histogram, SparseHistogram,
                                       HistogramPlotter2D)

class MockAxes(object):
    def __init__(self, xticks, yticks):
        self.xticks = xticks
        self.yticks = yticks

    def get_xticks(self): return self.xticks
    def get_yticks(self): return self.yticks

class TestHistogram(object):
    def setup(self):
        self.data = [1.0, 1.1, 1.2, 1.3, 2.0, 1.4, 2.3, 2.5, 3.1, 3.5]
        self.nbins = 5
        hist_counts = [5, 0, 2, 1, 1, 1]
        self.bins = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5]
        self.left_bin_edges = (1.0,)
        self.bin_widths = (0.5,)
        self.hist = collections.Counter({(0,): 5, (2,): 2, (3,): 1,
                                         (4,): 1, (5,): 1})

        self.default_hist = Histogram()
        self.hist_nbins = Histogram(n_bins=5)
        self.hist_nbins_binwidth = Histogram(n_bins=5, bin_width=1.0)
        self.hist_nbins_range = Histogram(n_bins=5, bin_range=(1.0, 3.5))
        self.hist_binwidth_range = Histogram(bin_width=0.5, bin_range=(1.0, 3.5))

    def test_initialization(self):
        assert_equal(self.default_hist.bins, 20)

        assert_equal(self.hist_nbins.bins, 5)
        assert_equal(self.hist_nbins.bin_width, None)

        assert_equal(self.hist_nbins_binwidth.bins, 5)
        assert_equal(self.hist_nbins_binwidth.bin_width, None)

        assert_items_equal(self.hist_nbins_range.bins, self.bins)
        assert_items_equal(self.hist_binwidth_range.bins, self.bins)


    def test_build_from_data(self):
        hist = self.hist_nbins.histogram(self.data)
        assert_equal(self.hist_nbins.count, 10)
        assert_items_equal(hist, self.hist)

        hist2 = self.hist_binwidth_range.histogram(self.data)
        assert_equal(self.hist_binwidth_range.count, 10)
        assert_items_equal(hist2, self.hist)

    @raises(RuntimeError)
    def test_build_from_data_fail(self):
        histo = Histogram(n_bins=5)
        histo.histogram()

    def test_add_data_to_histogram(self):
        histogram = Histogram(n_bins=5, bin_range=(1.0, 3.5))
        hist = histogram.add_data_to_histogram(self.data)
        assert_equal(histogram.count, 10)
        assert_items_equal(hist, self.hist)

        hist2 = histogram.add_data_to_histogram(self.data)
        assert_items_equal(hist2, hist+hist)
        assert_equal(histogram.count, 20)

    def test_compare_parameters(self):
        assert_equal(self.hist_nbins.compare_parameters(None), False)
        assert_equal(
            self.hist_nbins_range.compare_parameters(self.hist_binwidth_range),
            True
        )
        assert_equal(
            self.hist_binwidth_range.compare_parameters(self.hist_nbins_range),
            True
        )
        histo = Histogram(n_bins=5)
        assert_equal(self.hist_nbins_range.compare_parameters(histo), False)
        histo.histogram(self.data)
        assert_equal(self.hist_nbins_range.compare_parameters(histo), False)
        assert_equal(
            self.hist_nbins_range.compare_parameters(self.hist_nbins),
            False
        )
        assert_equal(histo.compare_parameters(self.hist_nbins), False)
        assert_equal(self.hist_nbins.compare_parameters(histo), False)

    def test_xvals(self):
        histo = Histogram(n_bins=5)
        hist = histo.histogram(self.data) # need this to set the bins
        assert_equal(histo.left_bin_edges, self.left_bin_edges)
        assert_equal(histo.bin_widths, self.bin_widths)
        assert_items_equal(histo.xvals("l"), [1.0, 1.5, 2.0, 2.5, 3.0, 3.5])
        assert_items_equal(histo.xvals("r"), [1.5, 2.0, 2.5, 3.0, 3.5, 4.0])
        assert_items_equal(histo.xvals("m"),
                           [1.25, 1.75, 2.25, 2.75, 3.25, 3.75])


    def test_normalization(self):
        histo = Histogram(n_bins=5)
        hist = histo.histogram(self.data)
        assert_equal(histo._normalization(), 5.0)

    def test_normalized(self):
        histo = Histogram(n_bins=5)
        hist = histo.histogram(self.data)
        assert_items_equal(list(histo.normalized().values()),
                           [1.0, 0.0, 0.4, 0.2, 0.2, 0.2])
        assert_items_equal(list(histo.normalized(raw_probability=True).values()),
                           [0.5, 0.0, 0.2, 0.1, 0.1, 0.1])

    def test_cumulative(self):
        histo = Histogram(n_bins=5)
        hist = histo.histogram(self.data)
        cumulative = list(histo.cumulative(None).values())
        assert_items_almost_equal(cumulative, [5.0, 5.0, 7.0, 8.0, 9.0, 10.0])
        assert_items_almost_equal(histo.cumulative(maximum=1.0),
                                  [0.5, 0.5, 0.7, 0.8, 0.9, 1.0])

    def test_cumulative_all_zero_warn(self):
        histo = Histogram(bin_width=0.5, bin_range=(1.0, 3.5))
        histo._histogram = collections.Counter({(0,): 0, (1,): 0})
        with pytest.warns(UserWarning, match=r"No non-zero"):
            histo.cumulative()

    def test_reverse_cumulative(self):
        histo = Histogram(n_bins=5)
        histo.histogram(self.data)
        rev_cumulative = histo.reverse_cumulative(maximum=None)
        assert_items_almost_equal(list(rev_cumulative.values()),
                                  [10, 5, 5, 3, 2, 1])
        rev_cumulative = histo.reverse_cumulative(maximum=1.0)
        assert_items_almost_equal(list(rev_cumulative.values()),
                                  [1.0, 0.5, 0.5, 0.3, 0.2, 0.1])

    def test_reverse_cumulative_all_zero_warn(self):
        histo = Histogram(bin_width=0.5, bin_range=(1.0, 3.5))
        histo._histogram = collections.Counter({(0,): 0, (1,): 0})
        with pytest.warns(UserWarning, match=r"No non-zero"):
            histo.reverse_cumulative()

    def test_left_bin_error(self):
        histo = Histogram(bin_width=0.5, bin_range=(-1.0, 3.5))
        histo.histogram([3.5])
        assert histo.reverse_cumulative() != 0

class TestSparseHistogram(object):
    def setup(self):
        data = [(0.0, 0.1), (0.2, 0.7), (0.3, 0.6), (0.6, 0.9)]
        self.histo = SparseHistogram(bin_widths=(0.5, 0.3),
                                     left_bin_edges=(0.0, -0.1))
        self.histo.histogram(data)

    def test_correct(self):
        correct_results = collections.Counter({
            (0, 0) : 1,
            (0, 2) : 2,
            (1, 3) : 1
        })
        assert_equal(self.histo._histogram, correct_results)

    def test_call(self):
        histo_fcn = self.histo()
        # voxels we have filled
        assert_equal(histo_fcn((0.25, 0.65)), 2)
        assert_equal(histo_fcn((0.01, 0.09)), 1)
        assert_equal(histo_fcn((0.61, 0.89)), 1)
        # empty voxel gives 0
        assert_equal(histo_fcn((2.00, 2.00)), 0)

    def test_normalized(self):
        raw_prob_normed = self.histo.normalized(raw_probability=True)
        assert_almost_equal(raw_prob_normed((0.25, 0.65)), 0.5)
        assert_almost_equal(raw_prob_normed((0.01, 0.09)), 0.25)
        assert_almost_equal(raw_prob_normed((0.61, 0.89)), 0.25)
        normed_fcn = self.histo.normalized()
        assert_almost_equal(normed_fcn((0.25, 0.65)), old_div(0.5,0.15))
        assert_almost_equal(normed_fcn((0.01, 0.09)), old_div(0.25,0.15))
        assert_almost_equal(normed_fcn((0.61, 0.89)), old_div(0.25,0.15))


class TestHistogramPlotter2D(object):
    def setup(self):
        data = [(0.0, 0.1), (0.2, 0.7), (0.3, 0.6), (0.6, 0.9)]
        histo = SparseHistogram(bin_widths=(0.5, 0.3),
                                left_bin_edges=(0.0, -0.1))
        histo.histogram(data)
        self.plotter = HistogramPlotter2D(histo)

    def test_to_bins(self):
        vals = [-0.1, 0.5, 0.8]
        xbins = self.plotter.to_bins(vals, 0)
        assert_items_almost_equal(xbins, [-0.2, 1.0, 1.6])
        ybins = self.plotter.to_bins(vals, 1)
        assert_items_almost_equal(ybins, [0.0, 2.0, 3.0])
        assert_equal(self.plotter.to_bins(None, 0), None)

    def test_axis_input(self):
        xt, xr, xl = self.plotter.axis_input(hist=[-1.0, 1.0, 2.0],
                                             ticklabels=None,
                                             lims=None,
                                             dof=0)
        assert_equal(xt, None)
        assert_items_equal(xr, (-1.0, 2.0))
        assert_items_equal(xl, (0, 3))

        xt, xr, xl = self.plotter.axis_input(hist=[-1.0, 1.0, 2.0],
                                             ticklabels=[-1.0, 0.0, 1.0],
                                             lims=None,
                                             dof=0)
        assert_items_equal(xt, [-2.0, 0.0, 2.0])
        assert_items_equal(xr, (-2.0, 2.0))
        assert_items_equal(xl, (0, 4.0))

        xt, xr, xl = self.plotter.axis_input(hist=[-1.0, 1.0, 2.0],
                                             ticklabels=[-1.0, 0.0, 1.0],
                                             lims=(-2.5, 0.0),
                                             dof=0)
        assert_items_equal(xt, [-2.0, 0.0, 2.0])
        assert_items_equal(xr, (-5.0, 2.0))
        assert_items_equal(xl, (0.0, 5.0))

    def test_ticks_and_labels(self):
        # mock axes, make sure they work as expected
        fake_ax = MockAxes([-1.0, 0.0, 1.0], [-6.0, 0.0, 6.0])
        assert_items_equal(fake_ax.get_xticks(), [-1.0, 0.0, 1.0])
        assert_items_equal(fake_ax.get_yticks(), [-6.0, 0.0, 6.0])

        old_format = self.plotter.label_format
        self.plotter.label_format = "{:4.2f}"

        xticks, xlabels = self.plotter.ticks_and_labels(
            ticks=[-2.0, -1.0, 0.0, 1.0, 2.0], ax=fake_ax, dof=0
        )
        assert_items_almost_equal(xticks, [-2.0, -1.0, 0.0, 1.0, 2.0])
        assert_items_equal(xlabels, ["-1.00", "-0.50", "0.00", "0.50", "1.00"])

        yticks, ylabels = self.plotter.ticks_and_labels(
            ticks=None, ax=fake_ax, dof=1
        )
        assert_items_almost_equal(yticks, [-6.0, 0.0, 6.0])
        assert_items_equal(ylabels, ["-1.90", "-0.10", "1.70"])

        self.plotter.label_format = old_format

