from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        raises, assert_almost_equal)
from nose.plugins.skip import SkipTest
from test_helpers import assert_items_almost_equal

import pandas as pd
import numpy as np
import openpathsampling as paths

import logging
logging.getLogger('openpathsampling.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.ensemble').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.netcdfplus').setLevel(logging.CRITICAL)


class testWHAM(object):
    def setup(self):
        self.exact = [1.0, 0.5, 0.25, 0.125, 0.0625, 0.3125]
        self.iface1 = [2.0, 1.0, 0.5, 0.25, 0.125, 0.0625, 0.0]
        self.iface2 = [1.0, 1.0, 1.0, 0.5, 0.25, 0.125, 0.0625]
        self.iface3 = [3.0, 3.0, 3.0, 3.0, 3.0, 1.5, 0.75]

        # self.iface1 = [1.0, 0.5, 0.25, 0.125, 0.0625, 0.0, 0.0]
        # self.iface2 = [1.0, 1.0, 1.0, 0.5, 0.25, 0.125, 0.0625]
        # self.iface3 = [1.0, 1.0, 1.0, 1.0, 1.0, 0.5, 0.25]

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
        self.wham = paths.analysis.WHAM(cutoff=0.1)


    def test_prep_reverse_cumulative(self):
        cleaned = self.wham.pandas_prep_reverse_cumulative(self.input_df)
        np.testing.assert_allclose(cleaned.as_matrix(),
                                   self.expected_cleaned)

    def test_unweighting_tis(self):
        unweighting = self.wham.pandas_unweighting_tis(self.cleaned)
        expected = np.array([[1.0, 0.0, 0.0],
                             [1.0, 0.0, 0.0],
                             [1.0, 1.0, 0.0],
                             [1.0, 1.0, 0.0],
                             [0.0, 1.0, 1.0],
                             [0.0, 1.0, 1.0],
                             [0.0, 0.0, 1.0]])
        np.testing.assert_allclose(unweighting.as_matrix(), expected)

    def test_sum_k_Hk_Q(self):
        sum_k_Hk_Q = self.wham.pandas_sum_k_Hk_Q(self.cleaned)
        expected = np.array([2.0, 1.0, 1.5, 0.75, 3.25, 1.625, 0.75])
        np.testing.assert_allclose(sum_k_Hk_Q.as_matrix(), expected)

    def test_n_entries(self):
        n_entries = self.wham.pandas_n_entries(self.cleaned)
        expected = np.array([3.75, 1.875, 5.25])
        np.testing.assert_allclose(n_entries.as_matrix(), expected)

    def test_weighted_counts_tis(self):
        n_entries = self.wham.pandas_n_entries(self.cleaned)
        unweighting = self.wham.pandas_unweighting_tis(self.cleaned)
        weighted_counts = self.wham.pandas_weighted_counts_tis(unweighting,
                                                               n_entries)
        expected = np.array([[3.75, 0.0, 0.0],
                             [3.75, 0.0, 0.0],
                             [3.75, 1.875, 0.0],
                             [3.75, 1.875, 0.0],
                             [0.0, 1.875, 5.25],
                             [0.0, 1.875, 5.25],
                             [0.0, 0.0, 5.25]])

        np.testing.assert_allclose(weighted_counts.as_matrix(), expected)

    def test_generate_lnZ(self):
        guess = [1.0, 1.0, 1.0]
        expected_lnZ = np.log([1.0, 1.0/4.0, 7.0/120.0])
        unweighting = self.wham.pandas_unweighting_tis(self.cleaned)
        sum_k_Hk_Q = self.wham.pandas_sum_k_Hk_Q(self.cleaned)
        weighted_counts = self.wham.pandas_weighted_counts_tis(
            unweighting,
            self.wham.pandas_n_entries(self.cleaned)
        )
        lnZ = self.wham.pandas_generate_lnZ(guess, unweighting,
                                            weighted_counts, sum_k_Hk_Q)
        np.testing.assert_allclose(lnZ.as_matrix(), expected_lnZ)
        raise SkipTest  # TODO: I'm not sure the last is log(7/120) 
        # (however, I got the same result out of the old version, too)



