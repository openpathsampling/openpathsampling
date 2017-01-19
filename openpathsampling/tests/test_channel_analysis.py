from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        raises, assert_almost_equal)
from nose.plugins.skip import SkipTest
from test_helpers import make_1d_traj

import openpathsampling as paths
import numpy as np
import random

import logging
logging.getLogger('openpathsampling.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.ensemble').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.netcdfplus').setLevel(logging.CRITICAL)

class testChannelAnalysis(object):
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
        self.toy_results =  {'a': [(0, 5), (8, 10)], 
                             'b': [(3, 9)],
                             'c': [(7, 9)]}
        self.set_a = frozenset(['a'])
        self.set_b = frozenset(['b'])
        self.set_c = frozenset(['c'])
        self.toy_expanded_results = [(0, 5, self.set_a), (3, 9, self.set_b), 
                                     (7, 9, self.set_c), (8, 10, self.set_a)]

    def _make_active(self, seq):
        traj = make_1d_traj(seq)
        sample = paths.Sample(replica=0,
                              trajectory=traj,
                              ensemble=self.ensemble)
        sample_set = paths.SampleSet(sample)
        return sample_set

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
        relabeled = paths.ChannelAnalysis._labels_by_step_newest(
            self.toy_expanded_results
        )
        assert_equal(relabeled, [(0, 3, self.set_a), (3, 7, self.set_b),
                                 (7, 8, self.set_c), (8, 10, self.set_a)])

    def test_labels_by_step_oldest(self):
        relabeled = paths.ChannelAnalysis._labels_by_step_oldest(
            self.toy_expanded_results
        )
        assert_equal(relabeled, [(0, 5, self.set_a), (5, 9, self.set_b),
                                 (9, 10, self.set_a)])

    def test_labels_by_step_multiple(self):
        relabeled = paths.ChannelAnalysis._labels_by_step_multiple(
            self.toy_expanded_results
        )
        assert_equal(relabeled, 
                     [(0, 3, self.set_a), 
                      (3, 5, self.set_a | self.set_b),
                      (5, 7, self.set_b),
                      (7, 8, self.set_b | self.set_c),
                      (8, 9, self.set_b | self.set_c | self.set_a),
                      (9, 10, self.set_a)])
