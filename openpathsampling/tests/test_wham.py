from __future__ import division
from __future__ import absolute_import
from builtins import object
from past.utils import old_div
from nose.tools import (assert_equal, assert_not_equal, raises,
                        assert_almost_equal)
from nose.plugins.skip import SkipTest
from .test_helpers import assert_items_almost_equal, assert_items_equal

import pandas as pd
import numpy as np
import openpathsampling as paths

import logging
logging.getLogger('openpathsampling.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.ensemble').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.netcdfplus').setLevel(logging.CRITICAL)


class TestWHAM(object):
    def setup(self):
        self.exact = [1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625]
        self.iface1 = [2.0, 1.0, 0.5, 0.25, 0.125, 0.0625, 0.0]
        self.iface2 = [1.0, 1.0, 1.0, 0.5, 0.25, 0.125, 0.0625]
        self.iface3 = [3.0, 3.0, 3.0, 3.0, 3.0, 1.5, 0.75]

        # self.iface1 = [1.0, 0.5, 0.25, 0.125, 0.0625, 0.0, 0.0]
        # self.iface2 = [1.0, 1.0, 1.0, 0.5, 0.25, 0.125, 0.0625]
        # self.iface3 = [1.0, 1.0, 1.0, 1.0, 1.0, 0.5, 0.25]

        # self.iface1 = [2.0, 0.5, 0.125, 0.0]
        # self.iface2 = [1.0, 1.0, 0.25, 0.0625]
        # self.iface3 = [3.0, 3.0, 3.0, 0.75]
        # self.index = [0.0, 0.2, 0.4, 0.6]

        self.columns = ["Interface 1", "Interface 2", "Interface 3"]
        self.index = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6]

        self.input_df = pd.DataFrame(
            data=np.array([self.iface1, self.iface2, self.iface3]).T,
            index=self.index,
            columns=self.columns
        )

        self.expected_cleaned = np.array([[2.0, 0.0, 0.0],
                                          [1.0, 0.0, 0.0],
                                          [0.5, 1.0, 0.0],
                                          [0.25, 0.5, 0.0],
                                          [0.0, 0.25, 3.0],
                                          [0.0, 0.125, 1.5],
                                          [0.0, 0.0, 0.75]])

        self.cleaned = pd.DataFrame(data=self.expected_cleaned,
                                    index=self.index,
                                    columns=self.columns)
        self.wham = paths.numerics.WHAM(cutoff=0.1)


    def test_prep_reverse_cumulative(self):
        cleaned = self.wham.prep_reverse_cumulative(self.input_df)
        np.testing.assert_allclose(cleaned.values,
                                   self.expected_cleaned)

    def test_prep_reverse_cumulative_with_interfaces(self):
        wham = paths.numerics.WHAM(cutoff=0.1, interfaces=[0.0, 0.2, 0.3])
        cleaned = wham.prep_reverse_cumulative(self.input_df)
        np.testing.assert_allclose(cleaned.values,
                                   np.array([[2.0, 0.0, 0.0],
                                             [1.0, 0.0, 0.0],
                                             [0.5, 1.0, 0.0],
                                             [0.25, 0.5, 3.0],
                                             [0.0, 0.25, 3.0],
                                             [0.0, 0.125, 1.5],
                                             [0.0, 0.0, 0.75]]))


    def test_unweighting_tis(self):
        unweighting = self.wham.unweighting_tis(self.cleaned)
        expected = np.array([[1.0, 0.0, 0.0],
                             [1.0, 0.0, 0.0],
                             [1.0, 1.0, 0.0],
                             [1.0, 1.0, 0.0],
                             [0.0, 1.0, 1.0],
                             [0.0, 1.0, 1.0],
                             [0.0, 0.0, 1.0]])
        np.testing.assert_allclose(unweighting.values, expected)

    def test_sum_k_Hk_Q(self):
        sum_k_Hk_Q = self.wham.sum_k_Hk_Q(self.cleaned)
        expected = np.array([2.0, 1.0, 1.5, 0.75, 3.25, 1.625, 0.75])
        np.testing.assert_allclose(sum_k_Hk_Q.values, expected)

    def test_n_entries(self):
        n_entries = self.wham.n_entries(self.cleaned)
        expected = np.array([3.75, 1.875, 5.25])
        np.testing.assert_allclose(n_entries.values, expected)

    def test_weighted_counts_tis(self):
        n_entries = self.wham.n_entries(self.cleaned)
        unweighting = self.wham.unweighting_tis(self.cleaned)
        weighted_counts = self.wham.weighted_counts_tis(unweighting,
                                                               n_entries)
        expected = np.array([[3.75, 0.0, 0.0],
                             [3.75, 0.0, 0.0],
                             [3.75, 1.875, 0.0],
                             [3.75, 1.875, 0.0],
                             [0.0, 1.875, 5.25],
                             [0.0, 1.875, 5.25],
                             [0.0, 0.0, 5.25]])

        np.testing.assert_allclose(weighted_counts.values, expected)

    def test_generate_lnZ(self):
        guess = [1.0, 1.0, 1.0]
        expected_lnZ = np.log([1.0, old_div(1.0,4.0), old_div(7.0,120.0)])
        # TODO: I'm not sure the last is log(7/120) 
        # however, I got the same result out of the old version, too, and
        # this does combine into the correct result in the end (see
        # test_output_histogram)
        unweighting = self.wham.unweighting_tis(self.cleaned)
        sum_k_Hk_Q = self.wham.sum_k_Hk_Q(self.cleaned)
        weighted_counts = self.wham.weighted_counts_tis(
            unweighting,
            self.wham.n_entries(self.cleaned)
        )
        lnZ = self.wham.generate_lnZ(guess, unweighting, weighted_counts,
                                     sum_k_Hk_Q)
        np.testing.assert_allclose(lnZ.values, expected_lnZ)

    def test_output_histogram(self):
        sum_k_Hk_Q = self.wham.sum_k_Hk_Q(self.cleaned)
        n_entries = self.wham.n_entries(self.cleaned)
        unweighting = self.wham.unweighting_tis(self.cleaned)
        weighted_counts = self.wham.weighted_counts_tis(unweighting,
                                                        n_entries)
        lnZ = pd.Series(data=np.log([1.0, old_div(1.0,4.0), old_div(7.0,120.0)]),
                        index=n_entries.index)
        wham_hist = self.wham.output_histogram(lnZ, sum_k_Hk_Q,
                                               weighted_counts)
        normed = self.wham.normalize_cumulative(wham_hist)
        np.testing.assert_allclose(normed.values, np.array(self.exact))


    def test_guess_lnZ_crossing_probability(self):
        input_data = np.array([[2.0, 1.0, 5.0],
                               [1.0, 1.0, 5.0],
                               [0.5, 1.0, 5.0],
                               [0.1, 0.2, 5.0],
                               [0.0, 0.04, 1.0],
                               [0.0, 0.02, 0.2]])
        input_df = pd.DataFrame(data=input_data,
                                index=self.index[0:6],
                                columns=self.columns)
        cleaned = self.wham.prep_reverse_cumulative(input_df)
        guess_lnZ = self.wham.guess_lnZ_crossing_probability(cleaned)
        expected_Z = np.array([1.0, 0.25, 0.25*0.2])
        np.testing.assert_allclose(guess_lnZ.values, np.log(expected_Z))

    def test_wham_bam_histogram(self):
        wham_hist = self.wham.wham_bam_histogram(self.input_df)
        np.testing.assert_allclose(wham_hist.values, self.exact)

    @raises(RuntimeError)
    def test_check_overlaps_no_overlap_with_first(self):
        bad_data = np.array([[1.0, 0.0, 0.0],
                             [0.5, 0.0, 0.0],
                             [0.0, 1.0, 0.0],
                             [0.0, 0.5, 1.0],
                             [0.0, 0.1, 0.2]])
        bad_df = pd.DataFrame(data=bad_data,
                              index=self.index[0:5],
                              columns=self.columns)
        self.wham.check_cleaned_overlaps(bad_df)

    @raises(RuntimeError)
    def test_check_overlaps_no_overlap_with_final(self):
        bad_data = np.array([[1.0, 0.0, 0.0],
                             [0.5, 0.0, 0.0],
                             [0.2, 1.0, 0.0],
                             [0.1, 0.5, 0.0],
                             [0.0, 0.0, 1.0],
                             [0.0, 0.0, 0.5]])
        bad_df = pd.DataFrame(data=bad_data,
                              index=self.index[0:6],
                              columns=self.columns)
        self.wham.check_cleaned_overlaps(bad_df)

    @raises(RuntimeError)
    def test_check_overlaps_no_overlap_in_middle(self):
        bad_data = np.array([[1.0, 0.0, 0.0, 0.0],
                             [0.5, 1.0, 0.0, 0.0],
                             [0.1, 0.2, 0.0, 0.0],
                             [0.0, 0.0, 1.0, 0.0],
                             [0.0, 0.0, 0.5, 1.0],
                             [0.0, 0.0, 0.1, 0.2]])
        bad_df = pd.DataFrame(data=bad_data,
                              index=self.index[0:6],
                              columns=self.columns + ['Interface 4'])
        self.wham.check_cleaned_overlaps(bad_df)
