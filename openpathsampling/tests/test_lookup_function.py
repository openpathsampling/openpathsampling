from __future__ import absolute_import
from builtins import object
from nose.tools import (assert_equal, assert_not_equal, assert_almost_equal,
                        raises, assert_in)
from nose.plugins.skip import SkipTest
from numpy import isnan
from .test_helpers import assert_items_almost_equal, assert_items_equal
import collections

import logging

from openpathsampling.numerics import (
    LookupFunction, LookupFunctionGroup, VoxelLookupFunction
)

class TestLookupFunctionGroup(object):
    def setup(self):
        x1 = [0.0, 1.0, 2.0, 3.0, 4.0]
        y1 = [0.1, 1.3, 1.8, 3.3, 4.2]
        luf1 = LookupFunction(x1, y1)
        x2 = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
        y2 = [0.6, 0.8, 1.3, 1.9, 2.6, 2.8, 3.3, 3.9]
        luf2 = LookupFunction(x2, y2)
        x3 = [1.0, 2.0]
        y3 = [1.0, 2.0]
        luf3 = LookupFunction(x3, y3)
        y4 = [1.1, 1.9]
        luf4 = LookupFunction(x3, y4)
        y5 = [1.2, 2.2]
        luf5 = LookupFunction(x3, y5)
        self.group = LookupFunctionGroup([luf1, luf2, luf3, luf4, luf5])

    def test_all_x(self):
        self.group.use_x = "all"
        assert_items_equal(
            self.group.x, [0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
        )

    def test_shared_x(self):
        self.group.use_x = "shared"
        assert_items_equal(self.group.x, [1.0, 2.0])

    def test_mean_shared(self):
        self.group.use_x = "shared"
        assert_almost_equal(self.group.mean(1.0), 1.08)
        assert_almost_equal(self.group(1.0), 1.08)

    def test_mean_allx(self):
        self.group.use_x = "all"
        assert_almost_equal(self.group.mean(1.5), 1.51)
        assert_almost_equal(self.group(1.5), 1.51)

    def test_std_shared(self):
        self.group.use_x = "shared"
        assert_almost_equal(self.group.std(1.0), 0.17204650534085253)

    def test_std_allx(self):
        self.group.use_x = "all"
        assert_almost_equal(self.group.std(1.5), 0.12806248474865695)

    def test_getitem(self):
        luf = self.group[0]
        assert_items_equal(luf.x, [0.0, 1.0, 2.0, 3.0, 4.0])
        assert_items_equal(luf, [0.1, 1.3, 1.8, 3.3, 4.2])

    def test_setitem(self):
        luf = LookupFunction([0.0, 1.0, 2.0], [0.5, 1.5, 2.5])
        self.group[3] = luf
        assert_almost_equal(self.group(1.0), 1.16)

    def test_append(self):
        luf = LookupFunction([0.0, 1.0, 2.0], [0.5, 1.5, 2.5])
        self.group.append(luf)
        assert_almost_equal(self.group(1.0), 1.15)


class TestVoxelLookupFunction(object):
    def setup(self):
        counter = collections.Counter({(0,0): 1.0, (0,2): 2.0, (1,4): 4.0,
                                       (-1,-1): 5.0})
        self.lookup = VoxelLookupFunction(left_bin_edges=(-1.0, 1.0),
                                          bin_widths=(0.5, 0.25),
                                          counter=counter)

    def test_keys(self):
        keys = [(0,0), (0,2), (1,4), (-1,-1)]
        assert_equal(len(list(self.lookup.keys())), len(keys))
        for key in list(self.lookup.keys()):
            assert_in(key, keys)

    def test_values(self):
        values = [1.0, 2.0, 4.0, 5.0]
        assert_equal(len(list(self.lookup.values())), len(values))
        for val in list(self.lookup.values()):
            assert_in(val, values)

    def test_bin_to_left_edge(self):
        assert_items_equal(self.lookup.bin_to_left_edge((0,0)), (-1.0, 1.0))
        assert_items_equal(self.lookup.bin_to_left_edge((2,-1)), (0.0, 0.75))

    def test_val_to_bin(self):
        result_bin = self.lookup.val_to_bin((-0.25, 2.0))
        assert_equal(len(result_bin), 2)
        assert_items_equal(result_bin, [1.5, 4.0])

    def test_counter_by_bin_edges(self):
        bin_edge_counter = self.lookup.counter_by_bin_edges
        expected = {(-1.0, 1.0): 1.0, (-1.0,1.5): 2.0, (-0.5, 2.0): 4.0,
                    (-1.5, 0.75): 5.0}
        assert_equal(len(expected), len(bin_edge_counter))
        for k in list(expected.keys()):
            assert_equal(bin_edge_counter[k], expected[k])

    def test_df_2d(self):
        df1 = self.lookup.df_2d()
        assert_items_equal(df1.index, [-1, 0, 1])
        assert_items_equal(df1.columns, [-1, 0, 2, 4])
        assert_equal(df1.at[-1, -1], 5.0)
        assert_equal(isnan(df1.at[0, 4]), True)

        df2 = self.lookup.df_2d(x_range=(-1, 3), y_range=(-1, 4))
        assert_items_equal(df2.index, [-1, 0, 1, 2, 3])
        assert_items_equal(df2.columns, [-1, 0, 1, 2, 3, 4])
        assert_equal(df2.at[-1, -1], 5.0)
        assert_equal(isnan(df2.at[0, 4]), True)
        assert_equal(isnan(df2.at[2, 3]), True)

    @raises(RuntimeError)
    def test_df_2d_not_2d(self):
        counter = collections.Counter({(0,0,0): 1.0})
        luf = VoxelLookupFunction(left_bin_edges=(0,0,0),
                                  bin_widths=(1,1,1),
                                  counter=counter)
        luf.df_2d()

    def test_call(self):
        assert_equal(self.lookup((0.0, 0.0)), 0.0)   # no bin
        assert_equal(self.lookup((-0.5, 2.0)), 4.0)  # bin edge
        assert_equal(self.lookup((-0.4, 2.1)), 4.0)  # in bin

