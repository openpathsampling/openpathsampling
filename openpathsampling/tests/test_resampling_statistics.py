import pandas as pd
from pandas.testing import assert_frame_equal
from nose.tools import assert_equal
from .test_helpers import assert_items_equal
import openpathsampling as paths

import logging
logging.getLogger('openpathsampling.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.ensemble').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.netcdfplus').setLevel(logging.CRITICAL)


class TestResamplingStatistics(object):
    # NOTE: we test the mean_df and std_df functions within this
    def setup(self):
        # order of the column/index labels should not matter
        self.list_AB = ['A', 'B']
        self.list_BA = ['B', 'A']
        self.df1 = pd.DataFrame([[1.5, 2.0], [2.1, 3.0]],
                                columns=['A', 'B'], index=['A', 'B'])
        self.df2 = pd.DataFrame([[2.5, 1.0], [3.5, 2.4]],
                                columns=['B', 'A'], index=['A', 'B'])
        self.df3 = pd.DataFrame([[2.25, 3.25], [1.25, 2.25]],
                                columns=['A', 'B'], index=['B', 'A'])
        self.df4 = pd.DataFrame([[3.25, 2.25], [2.25, 1.25]],
                                columns=['B', 'A'], index=['B', 'A'])
        self.inputs = [self.df1, self.df2, self.df3, self.df4]

    def test_std(self):
        # extra test to get the one line that isn't otherwise used
        from openpathsampling.numerics.resampling_statistics import std_df
        expected_std = pd.DataFrame(
            [[0.17677669529663689, 0.17677669529663689],
             [0.10606601717798207, 0.17677669529663689]],
            columns=['A', 'B'], index=['A', 'B']
        )
        assert_frame_equal(std_df(self.inputs), expected_std)

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
        expected_std = pd.DataFrame(
            [[0.17677669529663689, 0.17677669529663689],
             [0.10606601717798207, 0.17677669529663689]],
            columns=['A', 'B'], index=['A', 'B']
        )
        assert_frame_equal(stats.std, expected_std)

    def test_percentile(self):
        # TODO: this would benefit from more tests (with more input frames)
        stats = paths.numerics.ResamplingStatistics(
            function=lambda x: x,
            inputs=self.inputs
        )

        assert_frame_equal(
            stats.percentile(0),
            pd.DataFrame([[1.0, 2.0], [2.1, 3.0]],
                         index=self.list_AB, columns=self.list_AB),
            check_dtype=False  # better if this wasn't needed
        )

        assert_frame_equal(
            stats.percentile(100),
            pd.DataFrame([[1.5, 2.5], [2.4, 3.5]],
                         index=self.list_AB, columns=self.list_AB),
            check_dtype=False
        )

        assert_frame_equal(
            stats.percentile(50),
            pd.DataFrame([[1.25, 2.25], [2.25, 3.25]],
                         index=self.list_AB, columns=self.list_AB),
            check_dtype=False
        )


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
        resampler = paths.numerics.BlockResampling(self.samples, n_blocks=10)
        assert_equal(resampler.n_total_samples, 100)
        assert_equal(resampler.n_blocks, 10)
        assert_equal(resampler.n_per_block, 10)
        assert_equal(len(resampler.blocks), 10)
        assert_equal(resampler.n_resampled, 100)
        assert_equal(resampler.unassigned, [])
        assert_items_equal(resampler.blocks[0], list(range(0, 10)))
        assert_items_equal(resampler.blocks[1], list(range(10, 20)))
        assert_items_equal(resampler.blocks[5], list(range(50, 60)))
        assert_items_equal(resampler.blocks[-1], list(range(90, 100)))

        resampler_2 = paths.numerics.BlockResampling(self.samples, n_blocks=8)
        assert_equal(resampler_2.n_total_samples, 100)
        assert_equal(resampler_2.n_blocks, 8)
        assert_equal(resampler_2.n_per_block, 12)
        assert_equal(resampler_2.n_resampled, 96)
        assert_items_equal(resampler_2.blocks[-1], list(range(84, 96)))
        assert_items_equal(resampler_2.unassigned, list(range(96, 100)))

    def test_n_per_block_initialization(self):
        resampler = paths.numerics.BlockResampling(self.samples,
                                                   n_per_block=10)
        assert_equal(resampler.n_total_samples, 100)
        assert_equal(resampler.n_blocks, 10)
        assert_equal(resampler.n_per_block, 10)
        assert_equal(len(resampler.blocks), 10)
        assert_equal(resampler.n_resampled, 100)
        assert_equal(resampler.unassigned, [])
        assert_items_equal(resampler.blocks[0], list(range(0, 10)))
        assert_items_equal(resampler.blocks[1], list(range(10, 20)))
        assert_items_equal(resampler.blocks[5], list(range(50, 60)))
        assert_items_equal(resampler.blocks[-1], list(range(90, 100)))

        resampler_2 = paths.numerics.BlockResampling(self.samples,
                                                     n_per_block=8)
        assert_equal(resampler_2.n_total_samples, 100)
        assert_equal(resampler_2.n_blocks, 12)
        assert_equal(resampler_2.n_resampled, 96)
        assert_items_equal(resampler_2.blocks[-1], list(range(88, 96)))
        assert_items_equal(resampler_2.unassigned, list(range(96, 100)))

    def test_n_blocks_and_n_per_block_initialization(self):
        resampler = paths.numerics.BlockResampling(self.samples, n_blocks=9,
                                                   n_per_block=10)
        assert_equal(resampler.n_total_samples, 100)
        assert_equal(resampler.n_blocks, 9)
        assert_equal(resampler.n_per_block, 10)
        assert_equal(len(resampler.blocks), 9)
        assert_equal(resampler.n_resampled, 90)
        assert_items_equal(resampler.unassigned, list(range(90, 100)))
        assert_items_equal(resampler.blocks[0], list(range(0, 10)))
        assert_items_equal(resampler.blocks[1], list(range(10, 20)))
        assert_items_equal(resampler.blocks[5], list(range(50, 60)))
        assert_items_equal(resampler.blocks[-1], list(range(80, 90)))
