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

        columns = ["Interface 1", "Interface 2", "Interface 3"]
        index = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6]

        self.input_df = pd.DataFrame(
            data=np.array([self.iface1, self.iface2, self.iface3]).T,
            index=index,
            columns=columns
        )

        self.wham = paths.analysis.WHAM(cutoff=0.1)


    def test_prep_reverse_cumulative(self):
        cleaned = self.wham.pandas_prep_reverse_cumulative(self.input_df)
        expected = np.array([[2.0, 0.0, 0.0],
                             [1.0, 0.0, 0.0],
                             [0.5, 1.0, 0.0],
                             [0.25, 0.5, 0.0],
                             [0.0, 0.25, 3.0],
                             [0.0, 0.125, 1.5],
                             [0.0, 0.0, 0.75]])
        np.testing.assert_allclose(cleaned.as_matrix(), expected)

    def test_unweighting_tis(self):
        raise SkipTest

    def test_sum_k_Hk_Q(self):
        raise SkipTest

    def test_weighted_counts_tis(self):
        raise SkipTest

    def test_generate_lnZ(self):
        raise SkipTest

