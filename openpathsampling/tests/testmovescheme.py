from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        assert_almost_equal, assert_in, raises)
from nose.plugins.skip import Skip, SkipTest
from test_helpers import true_func, assert_equal_array_array, make_1d_traj

import openpathsampling as paths
from openpathsampling import VolumeFactory as vf
from openpathsampling.analysis.move_scheme import *
from openpathsampling.analysis.move_strategy import (
    MOVERLEVEL, GROUPLEVEL, SUPERGROUPLEVEL, GLOBALLEVEL,
    MoveStrategy, OneWayShootingStrategy, NearestNeighborRepExStrategy,
    DefaultStrategy
)

import logging
logging.getLogger('openpathsampling.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.ensemble').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)

class testMoveScheme(object):
    def setup(self):
        cvA = paths.CV_Function(name="xA", fcn=lambda s : s.xyz[0][0])
        cvB = paths.CV_Function(name="xA", fcn=lambda s : -s.xyz[0][0])
        self.stateA = paths.LambdaVolume(cvA, float("-inf"), -0.5)
        self.stateB = paths.LambdaVolume(cvB, float("-inf"), -0.5)
        interfacesA = vf.LambdaVolumeSet(cvA, float("-inf"), [-0.5, -0.3, 0.0])
        interfacesB = vf.LambdaVolumeSet(cvB, float("-inf"), [-0.5, -0.3, 0.0])
        network = paths.MSTISNetwork([
            (self.stateA, interfacesA, "A", cvA),
            (self.stateB, interfacesB, "B", cvB)
        ])
        self.scheme = MoveScheme(network)

    def test_sanity(self):
        raise SkipTest

    def test_append_individuals_default_levels(self):
        shootstrat = OneWayShootingStrategy()
        repexstrat = NearestNeighborRepExStrategy()
        defaultstrat = DefaultStrategy()
        assert_equal(len(self.scheme.strategies.keys()), 0)
        self.scheme.append(shootstrat)
        self.scheme.append(repexstrat)
        self.scheme.append(defaultstrat)

        strats = self.scheme.strategies
        assert_equal(len(strats.keys()), 3)
        for (k, v) in [(MOVERLEVEL, shootstrat), (GROUPLEVEL, repexstrat),
                       (GLOBALLEVEL, defaultstrat)]:
            assert_in(v, strats[k])

    def test_append_groups_default_levels(self):
        shootstrat = OneWayShootingStrategy()
        repexstrat = NearestNeighborRepExStrategy()
        defaultstrat = DefaultStrategy()
        assert_equal(len(self.scheme.strategies.keys()), 0)
        self.scheme.append([shootstrat, repexstrat, defaultstrat])

        strats = self.scheme.strategies
        assert_equal(len(strats.keys()), 3)
        for (k, v) in [(MOVERLEVEL, shootstrat), (GROUPLEVEL, repexstrat),
                       (GLOBALLEVEL, defaultstrat)]:
            assert_in(v, strats[k])


    def test_append_individuals_custom_levels(self):
        shootstrat = OneWayShootingStrategy()
        repexstrat = NearestNeighborRepExStrategy()
        defaultstrat = DefaultStrategy()
        assert_equal(len(self.scheme.strategies.keys()), 0)
        self.scheme.append(shootstrat, 60)
        self.scheme.append(repexstrat, 60)
        self.scheme.append(defaultstrat, 60)

        strats = self.scheme.strategies
        assert_equal(len(strats.keys()), 1)
        assert_items_equal(strats[60], [shootstrat, repexstrat, defaultstrat])

    def test_append_groups_same_custom_level(self):
        shootstrat = OneWayShootingStrategy()
        repexstrat = NearestNeighborRepExStrategy()
        defaultstrat = DefaultStrategy()
        assert_equal(len(self.scheme.strategies.keys()), 0)
        self.scheme.append([shootstrat, repexstrat, defaultstrat], 60)

        strats = self.scheme.strategies
        assert_equal(len(strats.keys()), 1)
        assert_items_equal(strats[60], [shootstrat, repexstrat, defaultstrat])

    def test_append_group_different_custom_levels(self):
        shootstrat = OneWayShootingStrategy()
        repexstrat = NearestNeighborRepExStrategy()
        defaultstrat = DefaultStrategy()
        assert_equal(len(self.scheme.strategies.keys()), 0)
        self.scheme.append([shootstrat, repexstrat, defaultstrat], 
                           [45, 55, 65])

        strats = self.scheme.strategies
        assert_equal(len(strats.keys()), 3)
        for (k, v) in [(45, shootstrat), (55, repexstrat), (65, defaultstrat)]:
            assert_in(v, strats[k])

    def test_include_movers(self):
        if self.scheme.movers == {}:
            print "Full support of MoveStrategy implemented?"
            print "Time to remove legacy from tests."
        else:
            self.scheme.movers = {} 

        shooters = self.scheme.network.movers['shooting']
        assert_equal(len(shooters), 4)

        self.scheme.include_movers(
            movers=shooters[:2], 
            group='shooting', 
            replace=False
        )
        raise SkipTest

    def test_default_move_decision_tree(self):
        raise SkipTest

