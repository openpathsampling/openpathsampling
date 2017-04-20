import pandas as pd
from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        assert_almost_equal, raises)
from nose.plugins.skip import Skip, SkipTest
from test_helpers import (
    true_func, assert_equal_array_array, make_1d_traj, data_filename
)

import openpathsampling as paths
from openpathsampling.high_level.interface_set import GenericVolumeInterfaceSet

import logging
logging.getLogger('openpathsampling.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.ensemble').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.netcdfplus').setLevel(logging.CRITICAL)

class TestResamplingStatistics(object):
    # NOTE: we test the mean_df and std_df functions within this
    def setup(self):
        pass

    def test_initialization(self):
        pass

    def test_percentile_range(self):
        pass

class TestBlockResampling(object):
    def setup(self):
        # create set of 100
        pass

    def test_default_initialization(self):
        pass

    def test_n_blocks_initialization(self):
        pass
