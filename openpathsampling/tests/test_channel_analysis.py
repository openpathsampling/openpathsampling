from __future__ import absolute_import
from builtins import object
from nose.tools import (assert_equal, assert_not_equal, raises,
                        assert_almost_equal)
from nose.plugins.skip import SkipTest
from numpy.testing import assert_array_almost_equal
from .test_helpers import make_1d_traj, data_filename, assert_items_equal

import openpathsampling as paths
import numpy as np
import random
import os

import logging
logging.getLogger('openpathsampling.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.ensemble').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.netcdfplus').setLevel(logging.CRITICAL)

class TestChannelAnalysis(object):
    def setup(self):
        cv = paths.FunctionCV("Id", lambda snap: snap.xyz[0][0])
        self.state_A = paths.CVDefinedVolume(cv, -0.1, 0.1)
        self.state_B = ~paths.CVDefinedVolume(cv, -1.0, 1.0)
        nml_increasing = paths.CVDefinedVolume(cv, 0.1, 1.0)
        nml_decreasing = paths.CVDefinedVolume(cv, -1.0, -0.1)
        increasing = paths.AllInXEnsemble(nml_increasing)
        decreasing = paths.AllInXEnsemble(nml_decreasing)
        self.ensemble = paths.SequentialEnsemble([
            paths.LengthEnsemble(1) & paths.AllInXEnsemble(self.state_A),
            paths.AllOutXEnsemble(self.state_A | self.state_B),
            paths.LengthEnsemble(1) & paths.AllInXEnsemble(self.state_B)
        ])
        self.incr_1 = self._make_active([0.0, 0.5, 1.1])
        self.incr_2 = self._make_active([0.05, 0.6, 1.2])
        self.decr_1 = self._make_active([0.0, -0.5, -1.1])
        self.both_1 = self._make_active([0.0, 0.5, -0.5, 1.1])
        self.both_2 = self._make_active([0.0, -0.4, 0.4, -1.1])
        self.none_1 = self._make_active([0.0, 1.1])
        self.none_2 = self._make_active([0.0, -1.1])

        self.channels = {
            'incr': increasing,
            'decr': decreasing
        }

        # used in simplest tests of relabeling
        self.toy_results = {'a': [(0, 5), (8, 10)],
                            'b': [(3, 9)],
                            'c': [(7, 9)]}
        self.results_with_none = {'a': [(0, 2), (6, 9)],
                                  'b': [(5, 7), (9, 10)],
                                  None: [(2, 5)]}
        self.set_a = frozenset(['a'])
        self.set_b = frozenset(['b'])
        self.set_c = frozenset(['c'])
        self.toy_expanded_results = [(0, 5, self.set_a), (3, 9, self.set_b),
                                     (7, 9, self.set_c), (8, 10, self.set_a)]
        self.expanded_results_simultaneous_ending = [
            (0, 5, self.set_a), (3, 9, self.set_b), (7, 10, self.set_c),
            (8, 10, self.set_a)
        ]
        self.expanded_oldest_skips_internal = [
            (0, 5, self.set_a), (3, 9, self.set_b), (7, 8, self.set_c),
            (8, 10, self.set_a), (10, 11, self.set_b)
        ]

    def _make_active(self, seq):
        traj = make_1d_traj(seq)
        sample = paths.Sample(replica=0,
                              trajectory=traj,
                              ensemble=self.ensemble)
        sample_set = paths.SampleSet(sample)
        return sample_set

    def test_storage(self):
        analyzer = paths.ChannelAnalysis(steps=None, channels=self.channels)
        analyzer._results = self.toy_results
        analyzer.treat_multiples = 'newest'
        storage = paths.Storage(data_filename('test.nc'), 'w')
        storage.save(self.incr_1[0])  # template snapshot
        storage.tag['analyzer'] = analyzer
        storage.sync()
        storage.close()

        new_store = paths.Storage(data_filename('test.nc'), 'r')
        reloaded = storage.tag['analyzer']
        assert_equal(reloaded._results, self.toy_results)

        if os.path.isfile(data_filename('test.nc')):
            os.remove(data_filename('test.nc'))

    def test_analyze_incr_decr(self):
        steps = [
            paths.MCStep(mccycle=0, active=self.incr_1),
            paths.MCStep(mccycle=1, active=self.decr_1)
        ]
        results = paths.ChannelAnalysis(steps, self.channels)
        assert_equal(results._results,
                     {'incr': [(0,1)], 'decr': [(1,2)], None: []})

    def test_analyze_incr_incr(self):
        steps = [
            paths.MCStep(mccycle=0, active=self.incr_1),
            paths.MCStep(mccycle=1, active=self.incr_2)
        ]
        results = paths.ChannelAnalysis(steps, self.channels)
        assert_equal(results._results,
                     {'incr': [(0,2)], 'decr': [], None: []})

    def test_analyze_incr_both_decr(self):
        steps = [
            paths.MCStep(mccycle=0, active=self.incr_1),
            paths.MCStep(mccycle=1, active=self.both_1),
            paths.MCStep(mccycle=2, active=self.both_2),
            paths.MCStep(mccycle=3, active=self.decr_1)
        ]
        results = paths.ChannelAnalysis(steps, self.channels)
        assert_equal(results._results,
                     {'incr': [(0,3)], 'decr': [(1,4)], None: []})

    def test_analyze_incr_both_incr(self):
        steps = [
            paths.MCStep(mccycle=0, active=self.incr_1),
            paths.MCStep(mccycle=1, active=self.both_1),
            paths.MCStep(mccycle=2, active=self.incr_2)
        ]
        results = paths.ChannelAnalysis(steps, self.channels)
        assert_equal(results._results,
                     {'incr': [(0,3)], 'decr': [(1,2)], None: []})

    def test_analyze_incr_same_decr(self):
        steps = [
            paths.MCStep(mccycle=0, active=self.incr_1),
            paths.MCStep(mccycle=1, active=self.incr_1),
            paths.MCStep(mccycle=2, active=self.incr_2)
        ]
        results = paths.ChannelAnalysis(steps, self.channels)
        assert_equal(results._results,
                     {'incr': [(0,3)], 'decr': [], None: []})

    def test_analyze_incr_none_incr(self):
        steps = [
            paths.MCStep(mccycle=0, active=self.incr_1),
            paths.MCStep(mccycle=1, active=self.none_1),
            paths.MCStep(mccycle=2, active=self.incr_2)
        ]
        results = paths.ChannelAnalysis(steps, self.channels)
        assert_equal(results._results,
                     {'incr': [(0,1), (2,3)], 'decr': [], None: [(1,2)]})

    def test_analyze_incr_none_decr(self):
        steps = [
            paths.MCStep(mccycle=0, active=self.incr_1),
            paths.MCStep(mccycle=1, active=self.none_1),
            paths.MCStep(mccycle=2, active=self.decr_1)
        ]
        results = paths.ChannelAnalysis(steps, self.channels)
        assert_equal(results._results,
                     {'incr': [(0,1)], 'decr': [(2,3)], None: [(1,2)]})

    def test_expand_results(self):
        expanded = paths.ChannelAnalysis._expand_results(self.toy_results)
        assert_equal(expanded, self.toy_expanded_results)

    def test_labels_by_step_newest(self):
        test_set = self.toy_expanded_results
        relabeled = paths.ChannelAnalysis._labels_by_step_newest(test_set)
        assert_equal(relabeled, [(0, 3, self.set_a), (3, 7, self.set_b),
                                 (7, 8, self.set_c), (8, 10, self.set_a)])

        test_set = self.expanded_results_simultaneous_ending
        relabeled = paths.ChannelAnalysis._labels_by_step_newest(test_set)
        assert_equal(relabeled, [(0, 3, self.set_a), (3, 7, self.set_b),
                                 (7, 8, self.set_c), (8, 10, self.set_a)])

    def test_labels_by_step_oldest(self):
        test_set = self.toy_expanded_results
        relabeled = paths.ChannelAnalysis._labels_by_step_oldest(test_set)
        assert_equal(relabeled, [(0, 5, self.set_a), (5, 9, self.set_b),
                                 (9, 10, self.set_a)])

        test_set = self.expanded_results_simultaneous_ending
        relabeled = paths.ChannelAnalysis._labels_by_step_oldest(test_set)
        assert_equal(relabeled, [(0, 5, self.set_a), (5, 9, self.set_b),
                                 (9, 10, self.set_c)])

        test_set = self.expanded_oldest_skips_internal
        relabeled = paths.ChannelAnalysis._labels_by_step_oldest(test_set)
        assert_equal(relabeled, [(0, 5, self.set_a), (5, 9, self.set_b),
                                 (9, 10, self.set_a), (10, 11, self.set_b)])


    def test_labels_by_step_multiple(self):
        test_set = self.toy_expanded_results
        relabeled = paths.ChannelAnalysis._labels_by_step_multiple(test_set)
        assert_equal(relabeled,
                     [(0, 3, self.set_a),
                      (3, 5, self.set_a | self.set_b),
                      (5, 7, self.set_b),
                      (7, 8, self.set_b | self.set_c),
                      (8, 9, self.set_b | self.set_c | self.set_a),
                      (9, 10, self.set_a)])

    def test_empty(self):
        obj = paths.ChannelAnalysis(steps=None, channels=self.channels)
        empty_results = {c: [] for c in list(self.channels.keys()) + [None]}
        assert_equal(obj._results, empty_results)

    def test_treat_multiples(self):
        obj = paths.ChannelAnalysis(steps=None, channels=self.channels)
        assert_equal(obj.treat_multiples, 'all')
        obj.treat_multiples = 'oldest'
        assert_equal(obj.treat_multiples, 'oldest')
        obj.treat_multiples = 'newest'
        assert_equal(obj.treat_multiples, 'newest')
        obj.treat_multiples = 'multiple'
        assert_equal(obj.treat_multiples, 'multiple')
        obj.treat_multiples = 'all'
        assert_equal(obj.treat_multiples, 'all')

    @raises(ValueError)
    def test_bad_treat_multiples(self):
        obj = paths.ChannelAnalysis(steps=None, channels=self.channels)
        assert_equal(obj.treat_multiples, 'all')
        obj.treat_multiples = 'bad_string'

    def test_labels_as_sets_sort_function(self):
        sort_key = paths.ChannelAnalysis._labels_as_sets_sort_function

        inp_list = [self.set_b, self.set_c, self.set_a]
        sorted_set_list = sorted(inp_list, key=sort_key)
        sorted_labels = [paths.ChannelAnalysis.label_to_string(e)
                         for e in sorted_set_list]
        assert_equal(sorted_labels, ['a', 'b', 'c'])

        inp_list = [self.set_a, self.set_a | self.set_b, self.set_b,
                    self.set_c | self.set_b,
                    self.set_a | self.set_b | self.set_c]
        sorted_set_list = sorted(inp_list, key=sort_key)
        sorted_labels = [paths.ChannelAnalysis.label_to_string(e)
                         for e in sorted_set_list]
        assert_equal(sorted_labels, ['a', 'b', 'a,b', 'b,c', 'a,b,c'])

    def test_labels_sort_function_with_none(self):
        sort_key = paths.ChannelAnalysis._labels_as_sets_sort_function
        inp_list = [self.set_b, self.set_a, set([None])]
        sorted_set_list = sorted(inp_list, key=sort_key)
        sorted_labels = [paths.ChannelAnalysis.label_to_string(e)
                         for e in sorted_set_list]
        assert_equal(sorted_labels, ['None', 'a', 'b'])

    def test_switching(self):
        analysis = paths.ChannelAnalysis(steps=None, channels=self.channels)
        analysis._results = self.toy_results
        #nan = float('nan')
        nan = 0  # self transitions are 0


        analysis.treat_multiples = 'newest'
        df = analysis.switching_matrix
        expected = np.array([[nan, 1, 0], [0, nan, 1], [1, 0, nan]])
        assert_array_almost_equal(df.values, expected)

        analysis.treat_multiples = 'oldest'
        df = analysis.switching_matrix
        expected = np.array([[nan, 1], [1, nan]]) # no column for c!
        assert_array_almost_equal(df.values, expected)

        analysis.treat_multiples = 'multiple'
        df = analysis.switching_matrix
        # 'a', 'b', 'a,b', 'b,c', 'a,b,c'
        expected = np.array([[nan, 0, 1, 0, 0],
                             [0, nan, 0, 1, 0],
                             [0, 1, nan, 0, 0],
                             [0, 0, 0, nan, 1],
                             [1, 0, 0, 0, nan]])
        assert_array_almost_equal(df.values, expected)

        # TODO: define switching when using 'all'

    def test_switching_with_none(self):
        analysis = paths.ChannelAnalysis(steps=None, channels=self.channels)
        analysis._results = self.results_with_none
        nan = 0  # self transitions are 0

        analysis.treat_multiples = 'newest'
        df = analysis.switching_matrix

        # None, a, b
        expected = np.array([[nan, 0, 1],
                             [1, nan, 1],
                             [0, 1, nan]])
        assert_array_almost_equal(df.values, expected)

    def test_residence_times(self):
        analysis = paths.ChannelAnalysis(steps=None, channels=self.channels)
        analysis._results = self.toy_results

        analysis.treat_multiples = 'newest'
        residence_times = analysis.residence_times
        assert_equal(residence_times, {'a': [3, 2], 'b': [4], 'c': [1]})

        analysis.treat_multiples = 'oldest'
        residence_times = analysis.residence_times
        assert_equal(residence_times, {'a': [5, 1], 'b': [4]})

        analysis.treat_multiples = 'all'
        residence_times = analysis.residence_times
        assert_equal(residence_times, {'a': [5, 2], 'b': [6], 'c': [2]})

        analysis.treat_multiples = 'multiple'
        residence_times = analysis.residence_times
        assert_equal(residence_times, {'a': [3, 1], 'b': [2], 'a,b': [2],
                                       'b,c': [1], 'a,b,c': [1]})

    def test_residence_times_with_none(self):
        analysis = paths.ChannelAnalysis(steps=None, channels=self.channels)
        analysis._results = self.results_with_none

        residence_times = analysis.residence_times
        assert_equal(residence_times,
                     {'a': [2, 3], 'b': [2, 1], 'None': [3]})

    def test_total_time(self):
        analysis = paths.ChannelAnalysis(steps=None, channels=self.channels)
        analysis._results = self.toy_results

        analysis.treat_multiples = 'newest'
        total_time = analysis.total_time
        assert_equal(total_time, {'a': 5, 'b': 4, 'c': 1})

        analysis.treat_multiples = 'oldest'
        total_time = analysis.total_time
        assert_equal(total_time, {'a': 6, 'b': 4})
        assert_equal(total_time['c'], 0)

        analysis.treat_multiples = 'all'
        total_time = analysis.total_time
        assert_equal(total_time, {'a': 7, 'b': 6, 'c': 2})

        analysis.treat_multiples = 'multiple'
        total_time = analysis.total_time
        assert_equal(total_time, {'a': 4, 'b': 2, 'a,b': 2, 'b,c': 1,
                                  'a,b,c': 1})

    def test_total_time_with_none(self):
        analysis = paths.ChannelAnalysis(steps=None, channels=self.channels)
        analysis._results = self.results_with_none

        total_time = analysis.total_time
        assert_equal(total_time, {'a': 5, 'b': 3, 'None': 3})


    def test_status(self):
        analysis = paths.ChannelAnalysis(steps=None, channels=self.channels)
        analysis._results = self.toy_results

        analysis.treat_multiples = 'newest'
        assert_equal(analysis.status(4), 'b')
        assert_equal(analysis.status(8), 'a')

        analysis.treat_multiples = 'oldest'
        assert_equal(analysis.status(4), 'a')
        assert_equal(analysis.status(8), 'b')

        analysis.treat_multiples = 'multiple'
        assert_equal(analysis.status(4), 'a,b')
        assert_equal(analysis.status(8), 'a,b,c')

        analysis.treat_multiples = 'all'
        assert_equal(analysis.status(4), 'a,b')
        assert_equal(analysis.status(8), 'a,b,c')

    def test_status_with_none(self):
        analysis = paths.ChannelAnalysis(steps=None, channels=self.channels)
        analysis._results = self.results_with_none
        assert_equal(analysis.status(4), 'None')

    @raises(RuntimeError)
    def test_bad_status_number(self):
        analysis = paths.ChannelAnalysis(steps=None, channels=self.channels)
        analysis._results = self.toy_results

        analysis.status(10)
