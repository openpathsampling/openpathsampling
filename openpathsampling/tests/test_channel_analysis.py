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

        self.channels = {
            'incr': increasing,
            'decr': decreasing
        }

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
        assert_equal(results._results, {'incr': [(0,1)], 'decr': [(1,2)]})

    def test_analyze_incr_incr(self):
        steps = [
            paths.MCStep(mccycle=0, active=self.incr_1),
            paths.MCStep(mccycle=1, active=self.incr_2)
        ]
        results = paths.ChannelAnalysis(steps, self.channels)
        assert_equal(results._results, {'incr': [(0,2)], 'decr': []})

    def test_analyze_incr_both_decr(self):
        steps = [
            paths.MCStep(mccycle=0, active=self.incr_1),
            paths.MCStep(mccycle=1, active=self.both_1),
            paths.MCStep(mccycle=2, active=self.both_2),
            paths.MCStep(mccycle=3, active=self.decr_1)
        ]
        results = paths.ChannelAnalysis(steps, self.channels)
        assert_equal(results._results, {'incr': [(0,3)], 'decr': [(1,4)]})

    def test_analyze_incr_both_incr(self):
        raise SkipTest

    def test_analyze_incr_same_decr(self):
        raise SkipTest

    def test_analyze_incr_none_incr(self):
        raise SkipTest

    def test_analyze_incr_none_decr(self):
        raise SkipTest
