from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        assert_almost_equal, raises)
from nose.plugins.skip import Skip, SkipTest
from test_helpers import true_func, assert_equal_array_array, make_1d_traj

import openpathsampling as paths
from openpathsampling.analysis.move_scheme import MoveScheme
from openpathsampling.analysis.move_strategy import *
from openpathsampling import VolumeFactory as vf

import logging
logging.getLogger('openpathsampling.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.ensemble').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)

class MoveStrategyTestSetup(object):
    def setup(self):
        cvA = paths.CV_Function(name="xA", fcn=lambda s : s.xyz[0][0])
        cvB = paths.CV_Function(name="xA", fcn=lambda s : -s.xyz[0][0])
        self.stateA = paths.LambdaVolume(cvA, float("-inf"), -0.5)
        self.stateB = paths.LambdaVolume(cvB, float("-inf"), -0.5)
        interfacesA = vf.LambdaVolumeSet(cvA, float("-inf"), [-0.5, -0.3, 0.0])
        interfacesB = vf.LambdaVolumeSet(cvB, float("-inf"), [-0.5, -0.3, 0.0])
        self.network = paths.MSTISNetwork([
            (self.stateA, interfacesA, "A", cvA),
            (self.stateB, interfacesB, "B", cvB)
        ])


class testMoveStrategy(MoveStrategyTestSetup):
    def test_get_ensembles(self):
        self.strategy = MoveStrategy(ensembles=None, network=self.network,
                                     group="test", replace=True)
        # load up the relevant ensembles to test against
        transition_ensembles = []
        for transition in self.network.sampling_transitions:
            transition_ensembles.append(transition.ensembles)
        assert_equal(len(transition_ensembles), 2)
        for ens_set in transition_ensembles:
            assert_equal(len(ens_set), 2)
        ensA = self.network.from_state[self.stateA].ensembles
        assert_equal(len(ensA), 2)
        # if you error before this, something is wrong in setup
        ensembles = self.strategy.get_ensembles(None)
        assert_equal(ensembles, transition_ensembles)

        ensembles = self.strategy.get_ensembles(ensA)
        assert_equal(ensembles, [ensA])

        extra_ens = transition_ensembles[1][0]
        weird_ens_list = [[ensA[0]], ensA[1], [extra_ens]]
        ensembles = self.strategy.get_ensembles(weird_ens_list)
        assert_equal(ensembles, [[ensA[0]], [ensA[1]], [extra_ens]])

class testOneWayShootingStrategy(MoveStrategyTestSetup):
    def test_make_movers(self):
        strategy = OneWayShootingStrategy()
        scheme = MoveScheme(self.network)
        movers = strategy.make_movers(scheme)
        assert_equal(len(movers), 4)
        for mover in movers:
            assert_equal(type(mover), paths.OneWayShootingMover)
            assert_equal(type(mover.selector), paths.UniformSelector)

class testNearestNeighborRepExStrategy(object):
    pass
