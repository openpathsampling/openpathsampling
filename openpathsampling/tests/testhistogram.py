from nose.tools import assert_equal, assert_not_equal, assert_items_equal, raises
from nose.plugins.skip import SkipTest
from test_helpers import assert_items_almost_equal

import logging
logging.getLogger('openpathsampling.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.ensemble').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.netcdfplus').setLevel(logging.CRITICAL)

import collections

from openpathsampling.analysis import Histogram

class testHistogram(object):
    def setup(self):
        self.data = [1.0, 1.1, 1.2, 1.3, 2.0, 1.4, 2.3, 2.5, 3.1, 3.5]
        self.nbins = 5
        hist_counts = [5, 0, 2, 1, 2]
        self.bins = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5]
        self.hist = collections.Counter({1.0: 5, 1.5: 0, 2.0: 2, 2.5: 1,
                                         3.0: 2})

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
        assert_equal(self.hist_nbins_range.compare_parameters(histo), True)
        assert_equal(
            self.hist_nbins_range.compare_parameters(self.hist_nbins),
            False
        )
        assert_equal(histo.compare_parameters(self.hist_nbins), True)
        assert_equal(self.hist_nbins.compare_parameters(histo), False)

    def test_xvals(self):
        histo = Histogram(n_bins=5)
        hist = histo.histogram(self.data) # need this to set the bins
        assert_items_equal(histo.bins, self.bins)
        assert_items_equal(histo.xvals("l"), [1.0, 1.5, 2.0, 2.5, 3.0])
        assert_items_equal(histo.xvals("r"), [1.5, 2.0, 2.5, 3.0, 3.5])
        assert_items_equal(histo.xvals("m"), [1.25, 1.75, 2.25, 2.75, 3.25])


    def test_normalization(self):
        histo = Histogram(n_bins=5)
        hist = histo.histogram(self.data)
        assert_equal(histo._normalization(), 5.0)

    def test_normalized(self):
        histo = Histogram(n_bins=5)
        hist = histo.histogram(self.data)
        assert_items_equal(histo.normalized(), [1.0, 0.0, 0.4, 0.2, 0.4])
        assert_items_equal(histo.normalized(raw_probability=True),
                           [0.5, 0.0, 0.2, 0.1, 0.2])

    def test_cumulative(self):
        histo = Histogram(n_bins=5)
        hist = histo.histogram(self.data)
        assert_items_almost_equal(histo.cumulative(), [5.0, 5.0, 7.0, 8.0, 10.0])
        assert_items_almost_equal(histo.cumulative(maximum=1.0), 
                                  [0.5, 0.5, 0.7, 0.8, 1.0])

    def test_reverse_cumulative(self):
        histo = Histogram(n_bins=5)
        hist = histo.histogram(self.data)
        assert_items_almost_equal(histo.reverse_cumulative(),
                                  [10, 5, 5, 3, 2])
        assert_items_almost_equal(histo.reverse_cumulative(maximum=1.0),
                                  [1.0, 0.5, 0.5, 0.3, 0.2])



