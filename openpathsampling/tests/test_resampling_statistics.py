import pandas as pd
from pandas.util.testing import assert_frame_equal
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
        # order of the column/index labels should not matter
        self.df1 = pd.DataFrame([[1.0, 2.0], [2.1, 3.0]],
                                columns=['A', 'B'], index=['A', 'B'])
        self.df2 = pd.DataFrame([[2.5, 1.5], [3.5, 2.4]],
                                columns=['B', 'A'], index=['A', 'B'])
        self.df3 = pd.DataFrame([[2.25, 3.25], [1.25, 2.25]],
                                columns=['A', 'B'], index=['B', 'A'])
        self.df4 = pd.DataFrame([[3.25, 2.25], [2.25, 1.25]],
                                columns=['B', 'A'], index=['B', 'A'])
        self.inputs = [self.df1, self.df2, self.df3, self.df4]

    def test_initialization(self):
        stats = paths.numerics.ResamplingStatistics(
            function=lambda x: x,
            inputs=self.inputs
        )

        for (truth, beauty) in zip(self.inputs, stats.results):
            assert_frame_equal(beauty, truth)

        expected_mean = pd.DataFrame([[1.25, 2.25], [2.25, 3.25]],
                                     columns=['A', 'B'], index=['A', 'B'])
        assert_frame_equal(stats.mean, expected_mean)
        # TODO: add test for std
        expected_std = pd.DataFrame(
            [[0.17677669529663689, 0.17677669529663689],
             [0.10606601717798207, 0.17677669529663689]],
            columns=['A', 'B'], index=['A', 'B']
        )
        assert_frame_equal(stats.std, expected_std)

    def test_percentile_range(self):
        raise SkipTest

class TestBlockResampling(object):
    def setup(self):
        self.samples = list(range(100))

    def test_default_initialization(self):
        resampler = paths.numerics.BlockResampling(self.samples)
        assert_equal(resampler.n_total_samples, 100)
        assert_equal(resampler.n_blocks, 20)
        assert_equal(resampler.n_per_block, 5)
        assert_equal(len(resampler.blocks), 20)
        assert_equal(resampler.n_resampled, 100)
        assert_equal(resampler.unassigned, [])
        assert_items_equal(resampler.blocks[0], list(range(0, 5)))
        assert_items_equal(resampler.blocks[1], list(range(5, 10)))
        assert_items_equal(resampler.blocks[10], list(range(50, 55)))
        assert_items_equal(resampler.blocks[-1], list(range(95, 100)))

    def test_n_blocks_initialization(self):
        raise SkipTest

    def test_n_per_block_initialization(self):
        raise SkipTest
