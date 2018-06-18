from __future__ import absolute_import
from builtins import object
from nose.tools import assert_equal, assert_not_equal, raises
from nose.plugins.skip import SkipTest
from .test_helpers import (
    CallIdentity, prepend_exception_message, make_1d_traj, data_filename,
    assert_items_equal
)


import openpathsampling as paths
from openpathsampling.high_level.transition import *

import logging
logging.getLogger('opentis.analysis.tis_analysis').setLevel(logging.DEBUG)
logging.getLogger('opentis.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.netcdfplus').setLevel(logging.CRITICAL)

class TestTISTransition(object):
    def setup(self):
        pass

    def test_initialization(self):
        pass

    def test_ensemble_statistics(self):
        pass

class TestFixedLengthTPSTransition(object):
    def setup(self):
        op = paths.FunctionCV("Id", lambda snap : snap.coordinates[0][0])
        self.stateA = paths.CVDefinedVolume(op, 0.1, 0.5)
        self.stateB = paths.CVDefinedVolume(op, 2.0, 2.5)

        self.transition = paths.FixedLengthTPSTransition(self.stateA,
                                                         self.stateB,
                                                         length=4)
        self.good_traj = self._make_traj("aixb")
        self.short_traj = self._make_traj("axb")
        self.long_traj = self._make_traj("aixbb")
        self.bad_states_traj = self._make_traj("aixx")

    def _make_traj(self, string):
        pretraj = []
        for l in string:
            pretraj.append({"a" : 0.3, "b" : 2.2, "i" : 0.6, "x" : 1.0}[l])
        return make_1d_traj(coordinates=pretraj, velocities=[1.0]*len(pretraj))

    def test_init(self):
        assert_equal(len(self.transition.ensembles), 1)
        assert_equal(self.transition.length, 4)
        assert_equal(self.transition.ensembles[0](self.good_traj), True)
        assert_equal(self.transition.ensembles[0](self.short_traj), False)
        assert_equal(self.transition.ensembles[0](self.long_traj), False)
        assert_equal(self.transition.ensembles[0](self.bad_states_traj), False)

    def test_storage(self):
        import os
        filename = data_filename("transitions.nc")
        if os.path.isfile(filename):
            os.remove(filename)
        storage = paths.Storage(filename, "w")
        storage.snapshots.save(self.good_traj[0])

        storage.save(self.transition)
        assert_equal(len(storage.transitions), 1)
        storage.sync_all()

        storage_r = paths.storage.AnalysisStorage(filename)
        reloaded = storage_r.transitions[0]
        assert_equal(self.transition.ensembles[0], reloaded.ensembles[0])

        if os.path.isfile(filename):
            os.remove(filename)

class TestMinusSidesSummary(object):
    def setup(self):
        op = paths.FunctionCV("Id", lambda snap : snap.coordinates[0][0])
        vol1 = paths.CVDefinedVolume(op, 0.1, 0.5)
        vol2 = paths.CVDefinedVolume(op, -0.1, 0.7)
        vol3 = paths.CVDefinedVolume(op, 2.0, 2.5)

        self.stateA = vol1
        self.innermost = vol2
        self.interstitial = vol2 & ~vol1
        self.outInterface = ~vol2 & ~vol3
        self.stateB = vol3

        # n_l = 2 tests
        self.traj_axaxa = self._make_traj("axxxxaaaaaxxxa")
        self.traj_aixiaixia = self._make_traj("aiixxxiaaaaaixiia")
        self.traj_aixixiaxia = self._make_traj("aiixxiixxxiaaaaaxxxiia")

        # n_l = 4 tests
        self.traj_axaxaxaxa = self._make_traj("axaaxxaaaxxxaaaaxxxxa")
        self.traj_aixixiaixaxiaixia = self._make_traj(
            "aixiixxiaaiixxaxiiaaaiixxxiiia"
        )

    def _make_traj(self, string):
        pretraj = []
        for l in string:
            pretraj.append({"a" : 0.3, "b" : 2.2, "i" : 0.6, "x" : 1.0}[l])
        return make_1d_traj(coordinates=pretraj, velocities=[1.0]*len(pretraj))

    def test_normal_minus(self):
        minus = paths.MinusInterfaceEnsemble(self.stateA, self.stateA)
        assert_equal(minus_sides_summary(self.traj_axaxa, minus),
                     {"in" : [5], "out" : [4]})
        assert_not_equal(minus_sides_summary(self.traj_axaxa, minus),
                         {"in" : [5], "out" : [3]})
        assert_equal(minus_sides_summary(self.traj_aixiaixia, minus), 
                     {"in" : [5], "out" : [6] })
        assert_equal(minus_sides_summary(self.traj_aixixiaxia, minus), 
                     {"in" : [5], "out" : [10] })

    def test_minus_with_interstitial(self):
        minus = paths.MinusInterfaceEnsemble(self.stateA, self.innermost)
        assert_equal(minus_sides_summary(self.traj_axaxa, minus),
                     {"in" : [5], "out" : [4]})
        assert_equal(minus_sides_summary(self.traj_aixiaixia, minus), 
                     {"in" : [6], "out" : [4] })
        assert_equal(minus_sides_summary(self.traj_aixixiaxia, minus), 
                     {"in" : [5], "out" : [8] })

    def test_minus_with_multiple_excursion(self):
        minus = paths.MinusInterfaceEnsemble(self.stateA, self.stateA,
                                             n_l=4)
        assert_equal(minus_sides_summary(self.traj_axaxaxaxa, minus),
                     {"in" : [2, 3, 4], "out" : [1, 2, 3]})
        assert_equal(minus_sides_summary(self.traj_aixixiaixaxiaixia, minus),
                     {"in" : [2, 1, 3], "out" : [7, 4, 3]})

    def test_minus_with_interstitial_and_multiple_excursion(self):
        minus = paths.MinusInterfaceEnsemble(self.stateA, self.innermost, n_l=4)
        assert_equal(minus_sides_summary(self.traj_axaxaxaxa, minus),
                     {"in" : [2, 3, 4], "out" : [1, 2, 3]})
        assert_equal(minus_sides_summary(self.traj_aixixiaixaxiaixia, minus),
                     {"in" : [4, 1, 5], "out" : [6, 2, 3]})
        assert_not_equal(minus_sides_summary(self.traj_aixixiaixaxiaixia, minus),
                     {"in" : [1, 4, 5], "out" : [6, 2, 3]})

