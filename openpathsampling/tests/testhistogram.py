from nose.tools import assert_equal, assert_not_equal, assert_items_equal, raises
from nose.plugins.skip import SkipTest

import logging

from openpathsampling.analysis import Histogram

class testHistogram(object):
    def setup(self):
        self.data = [1.0, 1.1, 1.2, 1.3, 2.0, 1.4, 2.3, 2.5, 3.1, 3.5]
        self.nbins = 5
        self.hist = [5, 0, 2, 1, 2]
        self.bins = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5]

        self.default_hist = Histogram()
        self.hist_nbins = Histogram(n_bins=5)
        self.hist_nbins_binwidth = Histogram(n_bins=5, bin_width=1.0)
        self.hist_nbins_range = Histogram(n_bins=5, bin_range=(1.0, 3.5))
        self.hist_binwidth_range = Histogram(bin_width=0.5, bin_range=(1.0, 3.5))

    def test_initialization(self):
        assert_equal(self.default_hist.bins, 40)

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

    def test_add_data_to_histogram(self):
        histogram = Histogram(n_bins=5, bin_range=(1.0, 3.5))
        hist = histogram.add_data_to_histogram(self.data)
        assert_equal(histogram.count, 10)
        assert_items_equal(hist, self.hist)

        hist2 = histogram.add_data_to_histogram(self.data)
        assert_items_equal(hist2, map((2).__mul__, hist))
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
        assert_equal(self.hist_nbins_range.compare_parameters(histo), True)
        assert_equal(
            self.hist_nbins_range.compare_parameters(self.hist_nbins),
            False
        )
        assert_equal(histo.compare_parameters(self.hist_nbins), True)
        assert_equal(self.hist_nbins.compare_parameters(histo), False)

    def test_normalized(self):
        raise SkipTest

    def test_cumulative(self):
        raise SkipTest

    def test_reverse_cumulative(self):
        raise SkipTest



