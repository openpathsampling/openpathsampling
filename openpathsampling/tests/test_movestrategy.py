from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from builtins import zip
from builtins import range
from past.utils import old_div
from builtins import object
from nose.tools import (assert_equal, assert_not_equal, assert_almost_equal,
                        raises, assert_in, assert_not_in)
from nose.plugins.skip import Skip, SkipTest
from .test_helpers import (
    true_func, assert_equal_array_array, make_1d_traj, MoverWithSignature,
    setify_ensemble_signature, reorder_ensemble_signature
)
import pytest

import openpathsampling as paths
from openpathsampling.high_level.move_scheme import MoveScheme, DefaultScheme
from openpathsampling.high_level.move_strategy import *

import collections

import logging
logging.getLogger('openpathsampling.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.ensemble').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.netcdfplus').setLevel(logging.CRITICAL)


class MockMoveStrategy(MoveStrategy):
    def make_movers(self, scheme):
        return None

class MockSingleEnsembleMoveStrategy(MockMoveStrategy,
                                     SingleEnsembleMoveStrategy):
    pass


def find_mover(scheme, group, sig):
    mover = None
    for m in scheme.movers[group]:
        if m.ensemble_signature_set == (set(sig[0]), set(sig[1])):
            mover = m
    return mover


class TestStrategyLevels(object):
    def test_level_type(self):
        assert_equal(levels.level_type(10), levels.SIGNATURE)
        assert_equal(levels.level_type(1), levels.SIGNATURE)
        assert_equal(levels.level_type(19), levels.SIGNATURE)
        assert_equal(levels.level_type(20), None)
        assert_equal(levels.level_type(21), levels.MOVER)
        assert_equal(levels.level_type(35), levels.MOVER)
        assert_equal(levels.level_type(100), levels.GLOBAL)


class MoveStrategyTestSetup(object):
    def setup(self):
        paths.InterfaceSet._reset()
        cvA = paths.FunctionCV(name="xA", f=lambda s: s.xyz[0][0])
        cvB = paths.FunctionCV(name="xB", f=lambda s: -s.xyz[0][0])
        self.stateA = paths.CVDefinedVolume(cvA, float("-inf"), -0.5)
        self.stateB = paths.CVDefinedVolume(cvB, float("-inf"), -0.5)
        interfacesA = paths.VolumeInterfaceSet(cvA, float("-inf"),
                                               [-0.5, -0.3, -0.1])
        interfacesB = paths.VolumeInterfaceSet(cvB, float("-inf"),
                                               [-0.5, -0.3, -0.1])
        self.network = paths.MSTISNetwork(
            [(self.stateA, interfacesA),
             (self.stateB, interfacesB)],
            ms_outers=paths.MSOuterTISInterface.from_lambdas(
                {interfacesA: 0.0, interfacesB: 0.0}
            )
        )


class TestMoveStrategy(MoveStrategyTestSetup):
    def test_levels(self):
        strategy = MockMoveStrategy(ensembles=None, group="test", replace=True)
        assert_equal(strategy.level, -1)
        assert_equal(strategy.replace_signatures, False)
        assert_equal(strategy.replace_movers, False)
        strategy.level = 10
        assert_equal(strategy.level, levels.SIGNATURE)
        assert_equal(strategy.replace_signatures, True)
        assert_equal(strategy.replace_movers, False)
        strategy.level = 25
        assert_not_equal(strategy.level, levels.MOVER)
        assert_equal(levels.level_type(strategy.level), levels.MOVER)
        assert_equal(strategy.replace_signatures, False)
        assert_equal(strategy.replace_movers, True)
        strategy.level = 99
        assert_equal(strategy.replace_signatures, False)
        assert_equal(strategy.replace_movers, False)

    def test_get_ensembles(self):
        self.strategy = MockMoveStrategy(ensembles=None, group="test",
                                         replace=True)
        scheme = MoveScheme(self.network)
        # load up the relevant ensembles to test against
        transition_ensembles = []
        for transition in self.network.sampling_transitions:
            transition_ensembles.append(transition.ensembles)
        assert_equal(len(transition_ensembles), 2)
        for ens_set in transition_ensembles:
            assert_equal(len(ens_set), 3)
        ensA = self.network.from_state[self.stateA].ensembles
        assert_equal(len(ensA), 3)
        # if you error before this, something is wrong in setup
        ensembles = self.strategy.get_ensembles(scheme, None)
        assert_equal(ensembles, transition_ensembles)

        ensembles = self.strategy.get_ensembles(scheme, ensA)
        assert_equal(ensembles, [ensA])

        extra_ens = transition_ensembles[1][0]
        weird_ens_list = [[ensA[0]], ensA[1], [extra_ens]]
        ensembles = self.strategy.get_ensembles(scheme, weird_ens_list)
        assert_equal(ensembles, [[ensA[0]], [ensA[1]], [extra_ens]])

        ensembles = self.strategy.get_ensembles(scheme, extra_ens)
        assert_equal(len(ensembles), 1)
        assert_equal(len(ensembles[0]), 1)
        assert_equal(ensembles[0][0], extra_ens)


class TestSingleEnsembleMoveStrategy(MoveStrategyTestSetup):
    def setup(self):
        super(TestSingleEnsembleMoveStrategy, self).setup()
        self.strategy = MockSingleEnsembleMoveStrategy(
            ensembles=None,
            group='test_group',
            replace=True
        )
        self.scheme = paths.DefaultScheme(self.network, engine=None)

    def test_get_per_mover_ensembles(self):
        per_mover_ensembles = \
                self.strategy.get_per_mover_ensembles(self.scheme)
        ensembles = self.scheme.network.sampling_ensembles
        assert len(per_mover_ensembles) == len(ensembles)
        for listed_ensemble in per_mover_ensembles:
            assert len(listed_ensemble) == 1
            assert listed_ensemble[0] in ensembles

    def test_get_parameters(self):
        list_params_as_list = [0, 1, 2, 3, 4, 5]
        list_params_as_single = 100
        nonlist_params = "alpha"
        params = self.strategy.get_parameters(
            scheme=self.scheme,
            list_parameters=[list_params_as_list, list_params_as_single],
            nonlist_parameters=[nonlist_params]
        )
        ensembles = self.scheme.network.sampling_ensembles
        expected = [(ens, num, list_params_as_single, nonlist_params)
                    for (ens, num) in zip(ensembles, list_params_as_list)]
        assert params == expected

    def test_get_parameters_error(self):
        with pytest.raises(RuntimeError):
            params = self.strategy.get_parameters(
                scheme=self.scheme,
                list_parameters=[[0, 1, 2, 3]]
            )


class TestForwardShootingStrategy(MoveStrategyTestSetup):
    def test_make_movers(self):
        strategy = ForwardShootingStrategy()
        scheme = MoveScheme(self.network)
        movers = strategy.make_movers(scheme)
        assert_equal(len(movers), 6)
        for mover in movers:
            assert_equal(type(mover), paths.ForwardShootMover)
            assert_equal(type(mover.selector), paths.UniformSelector)

    def test_make_movers_with_list(self):
        list_of_selectors = [
            paths.shooting.InterfaceConstrainedSelector(ens.interface)
            for ens in self.network.sampling_ensembles
        ]
        strategy = ForwardShootingStrategy(
            selector=list_of_selectors,
            ensembles=self.network.sampling_ensembles
        )
        scheme = MoveScheme(self.network)
        movers = strategy.make_movers(scheme)
        assert_equal(len(movers), 6)
        assert_equal(len(list_of_selectors), 6)
        for mover, sel in zip(movers, list_of_selectors):
            assert_equal(type(mover), paths.ForwardShootMover)
            assert_equal(type(mover.selector),
                         paths.shooting.InterfaceConstrainedSelector)
            assert_equal(mover.selector, sel)


class TestOneWayShootingStrategy(MoveStrategyTestSetup):
    def test_make_movers(self):
        strategy = OneWayShootingStrategy()
        scheme = MoveScheme(self.network)
        movers = strategy.make_movers(scheme)
        assert_equal(len(movers), 6)
        for mover in movers:
            assert_equal(type(mover), paths.OneWayShootingMover)
            assert_equal(type(mover.selector), paths.UniformSelector)

class TestTwoWayShootingStrategy(MoveStrategyTestSetup):
    def test_make_movers(self):
        strategy = TwoWayShootingStrategy(modifier=paths.NoModification())
        scheme = MoveScheme(self.network)
        movers = strategy.make_movers(scheme)
        assert_equal(len(movers), 6)
        for mover in movers:
            assert_equal(type(mover), paths.TwoWayShootingMover)
            assert_equal(type(mover.selector), paths.UniformSelector)
            assert_equal(type(mover.modifier), paths.NoModification)

    def test_composition_with_default_scheme(self):
        strategy = TwoWayShootingStrategy(modifier=paths.NoModification())
        scheme = DefaultScheme(self.network, engine=None)
        scheme.append(strategy)
        scheme.build_move_decision_tree()
        assert_equal(len(scheme.movers['shooting']), 6)
        for mover in scheme.movers['shooting']:
            assert_equal(type(mover), paths.TwoWayShootingMover)
            assert_equal(type(mover.selector), paths.UniformSelector)
            assert_equal(type(mover.modifier), paths.NoModification)


class TestNearestNeighborRepExStrategy(MoveStrategyTestSetup):
    def test_make_movers(self):
        strategy = NearestNeighborRepExStrategy()
        scheme = MoveScheme(self.network)
        movers = strategy.make_movers(scheme)
        assert_equal(len(movers), 4)
        ens0 = self.network.sampling_transitions[0].ensembles
        ens1 = self.network.sampling_transitions[1].ensembles
        assert_equal(movers[0].ensemble_signature_set,
                     (set([ens0[0], ens0[1]]), set([ens0[0], ens0[1]])))
        assert_equal(movers[1].ensemble_signature_set,
                     (set([ens0[1], ens0[2]]), set([ens0[1], ens0[2]])))
        assert_equal(movers[2].ensemble_signature_set,
                     (set([ens1[0], ens1[1]]), set([ens1[0], ens1[1]])))
        assert_equal(movers[3].ensemble_signature_set,
                     (set([ens1[1], ens1[2]]), set([ens1[1], ens1[2]])))

class TestAllSetRepExStrategy(MoveStrategyTestSetup):
    def test_make_movers(self):
        strategy = AllSetRepExStrategy()
        scheme = MoveScheme(self.network)
        movers = strategy.make_movers(scheme)
        assert_equal(len(movers), 6)
        ens0 = self.network.sampling_transitions[0].ensembles
        ens1 = self.network.sampling_transitions[1].ensembles

        signatures = [(set(m.ensemble_signature[0]),
                       set(m.ensemble_signature[1])) for m in movers]
        expected_signatures = [
            ((ens0[0], ens0[1]), (ens0[0], ens0[1])),
            ((ens0[0], ens0[2]), (ens0[0], ens0[2])),
            ((ens0[1], ens0[2]), (ens0[1], ens0[2])),
            ((ens1[0], ens1[1]), (ens1[0], ens1[1])),
            ((ens1[0], ens1[2]), (ens1[0], ens1[2])),
            ((ens1[1], ens1[2]), (ens1[1], ens1[2]))
        ]
        for sig in expected_signatures:
            set_sig = (set(sig[0]), set(sig[1]))
            assert_in(set_sig, signatures)

class TestSelectedPairsRepExStrategy(MoveStrategyTestSetup):
    def test_make_movers(self):
        ens00 = self.network.sampling_transitions[0].ensembles[0]
        ens02 = self.network.sampling_transitions[0].ensembles[2]
        strategy = SelectedPairsRepExStrategy(ensembles=[ens00, ens02])
        scheme = MoveScheme(self.network)
        movers = strategy.make_movers(scheme)
        assert_equal(len(movers), 1)
        assert_equal(movers[0].ensemble_signature_set,
                     ({ ens00, ens02 }, ({ ens00, ens02 })))

    @raises(RuntimeError)
    def test_init_ensembles_none(self):
        strategy = SelectedPairsRepExStrategy()

    @raises(RuntimeError)
    def test_init_ensembles_triplet(self):
        ensembles = self.network.sampling_transitions[0].ensembles
        strategy = SelectedPairsRepExStrategy(ensembles=ensembles)

    def test_make_movers_multiple_pairs(self):
        ens00 = self.network.sampling_transitions[0].ensembles[0]
        ens01 = self.network.sampling_transitions[0].ensembles[1]
        ens02 = self.network.sampling_transitions[0].ensembles[2]
        strategy = SelectedPairsRepExStrategy(ensembles=[[ens00, ens01],
                                                         [ens00, ens02],
                                                         [ens01, ens02]])
        scheme = MoveScheme(self.network)
        movers = strategy.make_movers(scheme)
        assert_equal(len(movers), 3)
        assert_equal(movers[0].ensemble_signature_set,
                     ({ens00, ens01}, {ens00, ens01}))
        assert_equal(movers[1].ensemble_signature_set,
                     ({ens00, ens02}, {ens00, ens02}))
        assert_equal(movers[2].ensemble_signature_set,
                     ({ens01, ens02}, {ens01, ens02}))

class TestReplicaExchangeStrategy(MoveStrategyTestSetup):
    def test_make_movers(self):
        strategy = ReplicaExchangeStrategy()
        scheme = MoveScheme(self.network)
        scheme.apply_strategy(NearestNeighborRepExStrategy())
        movers = strategy.make_movers(scheme)
        assert_equal(len(movers), 4)

    def test_swap_to_hop_to_swap(self):
        scheme = DefaultScheme(self.network)
        root = scheme.move_decision_tree()
        assert_equal(len(scheme.movers['repex']), 6)
        old_movers = scheme.movers['repex']
        scheme.append(EnsembleHopStrategy(), force=True)
        root = scheme.move_decision_tree(rebuild=True)
        assert_equal(len(scheme.movers['repex']), 12)
        scheme.append(ReplicaExchangeStrategy(), force=True)
        root = scheme.move_decision_tree(rebuild=True)
        assert_equal(len(scheme.movers['repex']), 6)
        new_movers = scheme.movers['repex']
        assert_not_equal(old_movers, new_movers)
        old_sigs = [m.ensemble_signature_set for m in old_movers]
        new_sigs = [m.ensemble_signature_set for m in new_movers]
        for old in old_sigs:
            assert_in(old, new_sigs)

    @raises(RuntimeError)
    def test_detailed_balance_partners(self):
        scheme = DefaultScheme(self.network)
        scheme.append(EnsembleHopStrategy(), force=True)
        root = scheme.move_decision_tree()
        assert_equal(len(scheme.movers['repex']), 12)
        scheme.movers['repex'].pop()
        assert_equal(len(scheme.movers['repex']), 11)
        strategy = ReplicaExchangeStrategy()
        strategy.make_movers(scheme)


class TestEnsembleHopStrategy(MoveStrategyTestSetup):
    def test_make_movers(self):
        strategy = EnsembleHopStrategy()
        scheme = MoveScheme(self.network)
        scheme.apply_strategy(NearestNeighborRepExStrategy())
        movers = strategy.make_movers(scheme)
        # defaults to 4 repex movers, so
        assert_equal(len(movers), 8)

        # set up the swap pairs
        swap_pairs = []
        for trans in self.network.sampling_transitions:
            swap_pairs.extend([[trans.ensembles[i], trans.ensembles[i+1]] 
                               for i in range(len(trans.ensembles)-1)])
        assert_equal(len(swap_pairs), 4)

        # check that each swap pair has a hop
        mover_sigs = [m.ensemble_signature for m in movers]
        for pair in swap_pairs:
            sig1 = ((pair[0],),(pair[1],))
            sig2 = ((pair[1],),(pair[0],))
            assert_in(sig1, mover_sigs)
            assert_in(sig2, mover_sigs)

        scheme.movers['repex'] = movers
        newmovers = strategy.make_movers(scheme)
        assert_equal(len(newmovers), 8)
        for mover in newmovers:
            assert_in(mover.ensemble_signature, mover_sigs)

    @raises(RuntimeError)
    def test_different_number_input_output_ensembles(self):
        ens0 = self.network.sampling_ensembles[0]
        ens1 = self.network.sampling_ensembles[1]
        ens2 = self.network.sampling_ensembles[2]
        weird_mover = MoverWithSignature(
            input_ensembles=[ens0, ens1, ens2],
            output_ensembles=[ens0, ens1]
        )
        assert_equal(weird_mover.ensemble_signature,
                     ((ens0,ens1,ens2),(ens0,ens1)))
        scheme = MoveScheme(self.network)
        scheme.movers['weird'] = [weird_mover]
        strategy = EnsembleHopStrategy(group='weird')
        strategy.make_movers(scheme)

    @raises(RuntimeError)
    def test_wrong_number_ensembles_in_signature(self):
        ens0 = self.network.sampling_ensembles[0]
        ens1 = self.network.sampling_ensembles[1]
        ens2 = self.network.sampling_ensembles[2]
        weird_mover = MoverWithSignature(
            input_ensembles=[ens0, ens1, ens2],
            output_ensembles=[ens0, ens1, ens2]
        )
        assert_equal(weird_mover.ensemble_signature,
                     ((ens0,ens1,ens2),(ens0,ens1,ens2)))
        scheme = MoveScheme(self.network)
        scheme.movers['weird'] = [weird_mover]
        strategy = EnsembleHopStrategy(group='weird')
        strategy.make_movers(scheme)

    @raises(RuntimeError)
    def test_not_replica_exchange_signature(self):
        ens0 = self.network.sampling_ensembles[0]
        ens1 = self.network.sampling_ensembles[1]
        ens2 = self.network.sampling_ensembles[2]
        weird_mover = MoverWithSignature(
            input_ensembles=[ens0, ens1],
            output_ensembles=[ens1, ens2]
        )
        assert_equal(weird_mover.ensemble_signature, 
                     ((ens0,ens1),(ens1,ens2)))
        scheme = MoveScheme(self.network)
        scheme.movers['weird'] = [weird_mover]
        strategy = EnsembleHopStrategy(group='weird')
        strategy.make_movers(scheme)

    def test_replace_nofrom(self):
        # this is the default behavior
        scheme = DefaultScheme(self.network)
        scheme.movers ={}
        scheme.append(EnsembleHopStrategy(replace=True, from_group=None))
        scheme.build_move_decision_tree()
        # 4 normal repex + 2 ms-outer repex = 6 repex * 2 hop/repex = 12
        assert_equal(len(scheme.movers['repex']), 12)

    def test_noreplace_from(self):
        # if replace is False and from_group is given, we end up with two
        # groups: both the new and the old.
        scheme = DefaultScheme(self.network)
        scheme.movers ={}
        scheme.append(EnsembleHopStrategy(replace=False, 
                                          group='hop',
                                          from_group='repex'))
        scheme.build_move_decision_tree()
        assert_equal(len(scheme.movers['repex']), 6)
        assert_equal(len(scheme.movers['hop']), 12)

    def test_replace_from(self):
        # if replace is True and we have a different from_group, we should
        # remove the old from_group from existence
        scheme = DefaultScheme(self.network)
        scheme.movers ={}
        scheme.append(EnsembleHopStrategy(replace=True, 
                                          group='hop',
                                          from_group='repex'))
        scheme.build_move_decision_tree()
        assert_equal(len(scheme.movers['hop']), 12)
        assert_not_in("repex", list(scheme.movers.keys()))

    def test_noreplace_nofrom(self):
        # if replace==False and the from_group is not given, we should
        # act as mover replacement (but this is seriously a bad idea)
        scheme = DefaultScheme(self.network)
        scheme.movers ={}
        scheme.append(EnsembleHopStrategy(replace=False, from_group=None))
        scheme.build_move_decision_tree()
        assert_equal(len(scheme.movers['repex']), 18)


class TestPathReversalStrategy(MoveStrategyTestSetup):
    def test_make_movers(self):
        strategy = PathReversalStrategy()
        scheme = MoveScheme(self.network)
        movers = strategy.make_movers(scheme)
        assert_equal(len(movers), 6)
        for m in movers:
            assert_equal(type(m), paths.PathReversalMover)


class TestMinusMoveStrategy(MoveStrategyTestSetup):
    def test_get_ensembles(self):
        strategy = MinusMoveStrategy()
        scheme = MoveScheme(self.network)
        ensembles = strategy.get_ensembles(scheme, None)
        assert_equal(len(ensembles), 2)
        for ens_group in ensembles:
            assert_equal(len(ens_group), 1)
        assert_not_equal(ensembles[0][0].state_vol, ensembles[1][0].state_vol)

    def test_get_ensembles_multiple_minus(self):
        strategy = MinusMoveStrategy()
        innerA = self.network.sampling_transitions[0].ensembles[0]
        innerB = self.network.sampling_transitions[1].ensembles[0]
        inner_vol_A = self.network.sampling_transitions[0].interfaces[0]
        inner_vol_B = self.network.sampling_transitions[1].interfaces[0]
        extra_minus = paths.MinusInterfaceEnsemble(
            state_vol=self.network.sampling_transitions[0].stateA,
            innermost_vols=[inner_vol_A, inner_vol_B]
        )
        self.network.special_ensembles['minus'][extra_minus] = [innerA, innerB]
        scheme = MoveScheme(self.network)
        ensembles = strategy.get_ensembles(scheme, None)
        assert_equal(len(ensembles), 2)
        assert_equal({ len(ensembles[0]), len(ensembles[1]) }, { 1, 2 })

    def test_get_ensembles_fixed_ensembles(self):
        strategy = MinusMoveStrategy()
        minusA = list(self.network.special_ensembles['minus'].keys())[0]
        scheme = MoveScheme(self.network)
        ensembles = strategy.get_ensembles(scheme, minusA)
        assert_equal(len(ensembles), 1)
        assert_equal(len(ensembles[0]), 1)
        assert_equal(ensembles[0][0], minusA)

    def test_make_movers(self):
        strategy = MinusMoveStrategy()
        scheme = MoveScheme(self.network)
        movers = strategy.make_movers(scheme)
        assert_equal(len(movers), 2)

        minuses = self.network.special_ensembles['minus']
        ens_minusA = list(minuses.keys())[0]
        ens_innerA = [t.ensembles[0] for t in minuses[ens_minusA]]
        sig_A = set([ens_minusA] + ens_innerA)
        ens_minusB = list(minuses.keys())[1]
        ens_innerB = [t.ensembles[0] for t in minuses[ens_minusB]]
        sig_B = set([ens_minusB] + ens_innerB)
        all_ens_sigs = [m.ensemble_signature_set for m in movers]

        # check the signatures
        assert_not_equal(sig_A, sig_B)
        assert_in(tuple([sig_A, sig_A]), all_ens_sigs)
        assert_in(tuple([sig_B, sig_B]), all_ens_sigs)

        # check that these are inner ensembles
        inners = [t.ensembles[0] for t in self.network.sampling_transitions]
        for inner in ens_innerA + ens_innerB:
            assert_in(inner, inners)

        # check that we've got the right inner for the right state
        stateA_inner = self.network.from_state[ens_minusA.state_vol].ensembles[0]
        assert_equal([stateA_inner], ens_innerA)
        stateB_inner = self.network.from_state[ens_minusB.state_vol].ensembles[0]
        assert_equal([stateB_inner], ens_innerB)

        # check that we've got minus ensembles
        for mover in movers:
            assert_in(mover.minus_ensemble, self.network.minus_ensembles)
            assert_equal(
                isinstance(mover.minus_ensemble, paths.MinusInterfaceEnsemble),
                True
            )


class TestSingleReplicaMinusMoveStrategy(MoveStrategyTestSetup):
    def test_make_movers(self):
        strategy = SingleReplicaMinusMoveStrategy()
        scheme = MoveScheme(self.network)
        movers = strategy.make_movers(scheme)
        assert_equal(len(movers), 2)

        minuses = self.network.special_ensembles['minus']
        ens_minusA = list(minuses.keys())[0]
        ens_innerA = [t.ensembles[0] for t in minuses[ens_minusA]]
        sig_A = set([ens_minusA] + ens_innerA)
        ens_minusB = list(minuses.keys())[1]
        ens_innerB = [t.ensembles[0] for t in minuses[ens_minusB]]
        sig_B = set([ens_minusB] + ens_innerB)
        all_ens_sigs = [m.ensemble_signature_set for m in movers]


class TestOrganizeByMoveGroupStrategy(MoveStrategyTestSetup):
    def scheme_setup_shooting_repex(self):
        scheme = MoveScheme(self.network)
        ens0 = self.network.sampling_transitions[0].ensembles[0]
        ens1 = self.network.sampling_transitions[0].ensembles[1]
        ens2 = self.network.sampling_transitions[0].ensembles[2]
        scheme.movers['shooting'] = [
            paths.OneWayShootingMover(
                selector=paths.UniformSelector(),
                ensemble=ens
            )
            for ens in [ens0, ens1, ens2]
        ]
        scheme.movers['repex'] = [
            paths.ReplicaExchangeMover(ens0, ens1),
            paths.ReplicaExchangeMover(ens1, ens2)
        ]
        return scheme

    def test_choice_probability(self):
        scheme = self.scheme_setup_shooting_repex()
        ens0 = self.network.sampling_transitions[0].ensembles[0]
        strategy = OrganizeByMoveGroupStrategy()
        group_weights = {'shooting' : 1.0, 'repex' : 0.5}
        mover_weights = {}
        for groupname in list(scheme.movers.keys()):
            for mover in scheme.movers[groupname]:
                mover_weights[(groupname, mover.ensemble_signature)] = 1.0

        choice_prob = strategy.choice_probability(scheme, group_weights,
                                                  mover_weights)
        # shooting moves: P_a*P_b = 1.0*1.0 (3 times)
        # repex moves: P_a*P_b = 0.5*1.0 (2 times)
        # norm = 3*1.0 + 2*0.5 = 4.0
        # each shooting prob: 0.25
        # each repex prob: 0.125
        assert_equal(len(choice_prob), len(sum(list(scheme.movers.values()), [])))
        for m in list(choice_prob.keys()):
            if m in scheme.movers['shooting']:
                assert_equal(choice_prob[m], 0.25)
            elif m in scheme.movers['repex']:
                assert_equal(choice_prob[m], 0.125)
            else:
                raise RuntimeError("Unknown mover "+repr(m))

        # now change the choice probability for one of them
        ens0_sig = ((ens0,),(ens0,))
        mover_weights[('shooting', ens0_sig)] = 2.0
        choice_prob = strategy.choice_probability(scheme, group_weights,
                                                  mover_weights)
        # shooting0: P_a*P_b = 1.0*2.0
        # new norm: 2*1.0 + 1*2.0 + 2*0.5 = 5.0
        # prob shooting0 : 2.0/5.0 = 0.4
        # prob shooting1,2 : 1.0/5.0 = 0.2
        # prob repex: 0.5/5.0 = 0.1
        for m in list(choice_prob.keys()):
            if m in scheme.movers['shooting']:
                   if m.ensemble_signature == ens0_sig:
                       assert_equal(choice_prob[m], 0.4)
                   else:
                       assert_equal(choice_prob[m], 0.2)
            elif m in scheme.movers['repex']:
                assert_equal(choice_prob[m], 0.1)
            else:
                raise RuntimeError("Unknown mover "+repr(m))

    @raises(KeyError)
    def test_choice_probability_bad_group_weights(self):
        scheme = MoveScheme(self.network)
        ens0 = self.network.sampling_transitions[0].ensembles[0]
        ens1 = self.network.sampling_transitions[0].ensembles[1]
        ens2 = self.network.sampling_transitions[0].ensembles[2]
        scheme.movers['shooting'] = [
            paths.OneWayShootingMover(
                selector=paths.UniformSelector(),
                ensemble=ens
            )
            for ens in [ens0, ens1, ens2]
        ]
        scheme.movers['repex'] = [
            paths.ReplicaExchangeMover(ens0, ens1),
            paths.ReplicaExchangeMover(ens1, ens2)
        ]
        strategy = OrganizeByMoveGroupStrategy()
        group_weights = {'shooting' : 1.0}
        mover_weights = {}
        for groupname in list(scheme.movers.keys()):
            for mover in scheme.movers[groupname]:
                mover_weights[(groupname, mover.ensemble_signature)] = 1.0
        assert_equal(len(mover_weights), 5)
        choice_prob = strategy.choice_probability(scheme, group_weights,
                                                  mover_weights)

    def test_weights_from_choice_probability(self):
        scheme = self.scheme_setup_shooting_repex()
        ens0 = self.network.sampling_transitions[0].ensembles[0]
        strategy = OrganizeByMoveGroupStrategy()
        group_weights = {'shooting' : 1.0, 'repex' : 0.5}
        mover_weights = {}
        for groupname in list(scheme.movers.keys()):
            for mover in scheme.movers[groupname]:
                mover_weights[(groupname, mover.ensemble_signature)] = 1.0
        ens0_sig = ((ens0,),(ens0,))
        mover_weights[('shooting', ens0_sig)] = 2.0

        choice_prob = strategy.choice_probability(scheme, group_weights,
                                                  mover_weights)

        (group_w, mover_w) = strategy.weights_from_choice_probability(
            scheme, choice_prob
        )

        assert_equal(group_weights, group_w)
        assert_equal(mover_weights, mover_w)
        #TODO: run more thorough tests of this

    def test_chooser_root_weights(self):
        scheme = MoveScheme(self.network)
        ens0 = self.network.sampling_transitions[0].ensembles[0]
        ens1 = self.network.sampling_transitions[0].ensembles[1]
        ens2 = self.network.sampling_transitions[0].ensembles[2]
        scheme.movers['shooting'] = [
            paths.OneWayShootingMover(
                selector=paths.UniformSelector(),
                ensemble=ens
            )
            for ens in [ens0, ens1, ens2]
        ]
        scheme.movers['repex'] = [
            paths.ReplicaExchangeMover(ens0, ens1),
            paths.ReplicaExchangeMover(ens1, ens2)
        ]
        strategy = OrganizeByMoveGroupStrategy()
        group_weights = {'shooting' : 1.0, 'repex' : 0.5}
        mover_weights = {}
        for groupname in list(scheme.movers.keys()):
            for mover in scheme.movers[groupname]:
                mover_weights[(groupname, mover.ensemble_signature)] = 1.0

        w = strategy.chooser_root_weights(scheme, group_weights, mover_weights)
        assert_equal(w, {'shooting' : 3.0, 'repex' : 1.0})

    def test_chooser_mover_weights(self):
        scheme = self.scheme_setup_shooting_repex()
        strategy = OrganizeByMoveGroupStrategy()
        group_weights = {'shooting' : 1.0, 'repex' : 0.5}
        mover_weights = {}
        for groupname in list(scheme.movers.keys()):
            for mover in scheme.movers[groupname]:
                mover_weights[(groupname, mover.ensemble_signature)] = 1.0
    
        m_sig = {}
        for g in scheme.movers:
            for m in scheme.movers[g]:
                m_sig[m] = (g, m.ensemble_signature)

        for g in scheme.movers:
            expected_w = {m : mover_weights[m_sig[m]] for m in scheme.movers[g]}
            w = strategy.chooser_mover_weights(scheme, g, mover_weights)
            assert_equal(w, expected_w)


    def test_make_movers(self):
        scheme = MoveScheme(self.network)
        ens0 = self.network.sampling_transitions[0].ensembles[0]
        ens1 = self.network.sampling_transitions[0].ensembles[1]
        ens2 = self.network.sampling_transitions[0].ensembles[2]
        scheme.movers['shooting'] = [
            paths.OneWayShootingMover(
                selector=paths.UniformSelector(),
                ensemble=ens
            )
            for ens in [ens0, ens1, ens2]
        ]
        scheme.movers['repex'] = [
            paths.ReplicaExchangeMover(ensemble1=ens0, ensemble2=ens1),
            paths.ReplicaExchangeMover(ensemble1=ens1, ensemble2=ens2)
        ]
        scheme.movers['pathreversal'] = [
            paths.PathReversalMover(ensemble=ens)
            for ens in [ens0, ens1, ens2]
        ]
        scheme.movers['minus'] = [paths.MinusMover(
            minus_ensemble=self.network.minus_ensembles[0],
            innermost_ensembles=[ens0]
        )]

        strategy = OrganizeByMoveGroupStrategy()
        root = strategy.make_movers(scheme)
        
        assert_equal(len(root.movers), 4)
        names = ['ShootingChooser', 'RepexChooser', 'PathreversalChooser', 
                 'MinusChooser']
        name_dict = {root.movers[i].name : i for i in range(len(root.movers))}
        for name in names:
            assert_in(name, list(name_dict.keys()))

        name = 'ShootingChooser'
        weight = root.weights[name_dict[name]]
        chooser = root.movers[name_dict[name]]
        assert_equal(type(chooser), paths.RandomChoiceMover)
        assert_equal(weight, 3.0)
        assert_equal(len(chooser.movers), 3)
        for w in chooser.weights:
            assert_equal(w, 1.0)

        name = 'RepexChooser'
        weight = root.weights[name_dict[name]]
        chooser = root.movers[name_dict[name]]
        assert_equal(type(chooser), paths.RandomChoiceMover)
        assert_equal(weight, 1.0)
        assert_equal(len(chooser.movers), 2)
        for w in chooser.weights:
            assert_equal(w, 1.0)

        name = 'MinusChooser'
        weight = root.weights[name_dict[name]]
        chooser = root.movers[name_dict[name]]
        assert_equal(type(chooser), paths.RandomChoiceMover)
        assert_equal(len(chooser.movers), 1)
        assert_equal(weight, 0.2)
        for w in chooser.weights:
            assert_equal(w, 1.0)

    def test_make_movers_unknown_group(self):
        scheme = MoveScheme(self.network)
        ens0 = self.network.sampling_transitions[0].ensembles[0]
        ens1 = self.network.sampling_transitions[0].ensembles[1]
        ens2 = self.network.sampling_transitions[0].ensembles[2]
        scheme.movers['blahblah']  = [
            paths.ReplicaExchangeMover(ensemble1=ens0, ensemble2=ens1),
            paths.ReplicaExchangeMover(ensemble1=ens1, ensemble2=ens2)
        ]

        strategy = OrganizeByMoveGroupStrategy()
        root = strategy.make_movers(scheme)

        name_dict = {root.movers[i].name : i for i in range(len(root.movers))}

        name = 'BlahblahChooser'
        weight = root.weights[name_dict[name]]
        chooser = root.movers[name_dict[name]]
        assert_equal(type(chooser), paths.RandomChoiceMover)
        assert_equal(weight, 2.0)
        assert_equal(len(chooser.movers), 2)
        for w in chooser.weights:
            assert_equal(w, 1.0)

    def test_make_movers_custom_group(self):
        scheme = MoveScheme(self.network)
        ens0 = self.network.sampling_transitions[0].ensembles[0]
        ens1 = self.network.sampling_transitions[0].ensembles[1]
        ens2 = self.network.sampling_transitions[0].ensembles[2]
        scheme.movers['blahblahblah']  = [
            paths.ReplicaExchangeMover(ensemble1=ens0, ensemble2=ens1),
            paths.ReplicaExchangeMover(ensemble1=ens1, ensemble2=ens2)
        ]

        strategy = OrganizeByMoveGroupStrategy()
        strategy.group_weights['blahblahblah'] = 2.0
        root = strategy.make_movers(scheme)

        name_dict = {root.movers[i].name : i for i in range(len(root.movers))}
        
        name = 'BlahblahblahChooser'
        weight = root.weights[name_dict[name]]
        chooser = root.movers[name_dict[name]]
        assert_equal(type(chooser), paths.RandomChoiceMover)
        assert_equal(weight, 4.0)
        assert_equal(len(chooser.movers), 2)
        for w in chooser.weights:
            assert_equal(w, 1.0)

    def test_get_weights_scheme_all_unset(self):
        strategy = OrganizeByMoveGroupStrategy()

        scheme = MoveScheme(self.network)
        scheme.append(NearestNeighborRepExStrategy())
        scheme.append(OneWayShootingStrategy())
        root = scheme.move_decision_tree()
        assert_equal(len(scheme.movers), 2)
        all_movers = scheme.movers['shooting'] + scheme.movers['repex']
        all_movers_sigs = [m.ensemble_signature for m in all_movers]
        assert_equal(strategy.group_weights, {})
        assert_equal(strategy.mover_weights, {})

        (group_weights, mover_weights) = strategy.get_weights(
            scheme=scheme,
            sorted_movers=scheme.movers
        )
        assert_equal(group_weights, {'shooting' : 1.0, 'repex' : 0.5})

        # check that the number of sigs in a each group matches the number
        # of movers in that group
        for group in list(group_weights.keys()):
            mover_sigs = [sig for sig in list(mover_weights.keys()) 
                          if sig[0]==group]
            assert_equal(len(mover_sigs), len(scheme.movers[group]))

        for sig in list(mover_weights.keys()):
            assert_equal(mover_weights[sig], 1.0)
            assert_in(sig[1], all_movers_sigs)

        # check that we can reuse these in a different scheme
        scheme2 = MoveScheme(self.network)
        scheme2.append(OneWayShootingStrategy())
        root = scheme2.move_decision_tree()
        assert_equal(len(scheme2.movers), 1)

        (group_weights, mover_weights) = strategy.get_weights(
            scheme=scheme2,
            sorted_movers=scheme2.movers
        )
        assert_equal(group_weights, {'shooting' : 1.0})
        for sig in mover_weights:
            assert_equal(mover_weights[sig], 1.0)
            assert_in(sig[1], [m.ensemble_signature 
                               for m in scheme2.movers[sig[0]]])

        for group in list(scheme2.movers.keys()):
            mover_sigs = [sig for sig in list(mover_weights.keys())
                          if sig[0]==group]
            assert_equal(len(mover_sigs), len(scheme2.movers[group]))

    def test_get_weights_both_internal_weights_set(self):
        strategy = OrganizeByMoveGroupStrategy()
        ensA = self.network.sampling_transitions[0].ensembles[0]
        ensA_sig = ((ensA,),(ensA,))

        scheme = MoveScheme(self.network)
        scheme.append([NearestNeighborRepExStrategy(), 
                       OneWayShootingStrategy(), 
                       strategy])
        root = scheme.move_decision_tree()
        old_choice_probability = scheme.choice_probability
        old_shooter_ensA = [m for m in scheme.movers['shooting'] 
                            if m.ensemble_signature == ensA_sig][0]

        strategy.group_weights['repex'] = 3.0
        strategy.mover_weights[('shooting',ensA_sig)]= 2.0
        scheme.choice_probability = {} # this is unset here
 
        (group_weights, mover_weights) = strategy.get_weights(
            scheme=scheme,
            sorted_movers=scheme.movers,
            sort_weights_override=strategy.group_weights,
            mover_weights_override=strategy.mover_weights
        )
        root = scheme.move_decision_tree(rebuild=True)

        assert_equal(group_weights, {'shooting' : 1.0, 'repex' : 3.0})
        expected_mover_weights = {}
        for group in scheme.movers:
            for mover in scheme.movers[group]:
                sig = mover.ensemble_signature
                if group == 'shooting' and sig == ensA_sig:
                    expected_mover_weights[(group,sig)] = 2.0
                elif group == 'shooting':
                    expected_mover_weights[(group,sig)] = 1.0
                else:
                    expected_mover_weights[(group,sig)] = 1.0

        assert_equal(mover_weights, expected_mover_weights)

        new_choice_probability = scheme.choice_probability
        new_repex_chooser = [m for m in root if m.name=="RepexChooser"][0]
        repex_chooser_idx = root.movers.index(new_repex_chooser)
        assert_equal(root.weights[repex_chooser_idx],
                     3.0*len(new_repex_chooser.movers))

        new_shoot_chooser = [m for m in root if m.name=="ShootingChooser"][0]
        new_shooter_ensA = [m for m in scheme.movers['shooting'] 
                            if m.ensemble_signature == ensA_sig][0]
        shooter_ensA_idx = new_shoot_chooser.movers.index(new_shooter_ensA)
        assert_equal(new_shoot_chooser.weights[shooter_ensA_idx], 2.0)

        assert_not_equal(new_choice_probability, old_choice_probability)

    def test_get_weights_group_weights_set(self):
        strategy = OrganizeByMoveGroupStrategy()
        scheme = MoveScheme(self.network)
        scheme.append([NearestNeighborRepExStrategy(), 
                       OneWayShootingStrategy(), 
                       strategy])
        root = scheme.move_decision_tree()

        nshoot = len(scheme.movers['shooting'])
        nrepex = len(scheme.movers['repex'])

        for mover in scheme.movers['shooting']:
            assert_almost_equal(scheme.choice_probability[mover], old_div(1.0,8.0))
        for mover in scheme.movers['repex']:
            assert_almost_equal(scheme.choice_probability[mover], old_div(1.0,16.0))


        strategy.group_weights['shooting'] = 2.0
        root = scheme.move_decision_tree(rebuild=True)

        assert_equal(strategy.group_weights, {'shooting' : 2.0, 'repex' : 0.5})

        (group_weights, mover_weights) = strategy.get_weights(
            scheme=scheme,
            sorted_movers=scheme.movers,
            sort_weights_override=strategy.group_weights
        )
        assert_equal(group_weights, {'shooting' : 2.0, 'repex' : 0.5})
        # everything within the group should have the same mover_weight
        for group in scheme.movers:
            group_sigs = [s for s in mover_weights if s[0]==group]
            for sig in group_sigs:
                assert_equal(mover_weights[sig], mover_weights[group_sigs[0]])

        choice_prob = strategy.choice_probability(
            scheme, group_weights, mover_weights
        )
        for mover in scheme.movers['shooting']:
            assert_almost_equal(choice_prob[mover], old_div(1.0,7.0))
        for mover in scheme.movers['repex']:
            assert_almost_equal(choice_prob[mover], old_div(1.0,28.0))
        

    def test_get_weights_mover_weights_set(self):
        strategy = OrganizeByMoveGroupStrategy()
        scheme = MoveScheme(self.network)
        scheme.append([NearestNeighborRepExStrategy(), 
                       OneWayShootingStrategy(),
                       strategy])
        # build it so there's an old version
        strategy.group_weights['repex'] = 3.0
        root = scheme.move_decision_tree()
        (group_weights, mover_weights) = strategy.get_weights(
            scheme=scheme,
            sorted_movers=scheme.movers,
            sort_weights_override=strategy.group_weights
        )
        assert_equal(group_weights, {'shooting' : 1.0, 'repex' : 3.0})

        ensA = self.network.sampling_transitions[0].ensembles[0]
        ensA_sig = ((ensA,),(ensA,)) 
        strategy.group_weights = {}
        strategy.mover_weights[('shooting',ensA_sig)] = 2.0

        (group_weights, mover_weights) = strategy.get_weights(
            scheme=scheme,
            sorted_movers=scheme.movers,
            sort_weights_override=strategy.group_weights,
            mover_weights_override=strategy.mover_weights
        )

        for sig in mover_weights:
            if sig == ('shooting',ensA_sig):
                assert_equal(mover_weights[sig], 2.0)
            elif sig[0] == 'shooting':
                assert_equal(mover_weights[sig], 1.0)

        assert_almost_equal(group_weights['shooting'], 1.0)
        assert_almost_equal(group_weights['repex'], 3.0)

    def test_get_weights_internal_unset_choice_prob_set(self):
        strategy = OrganizeByMoveGroupStrategy()
        scheme = MoveScheme(self.network)
        scheme.append([NearestNeighborRepExStrategy(), 
                       OneWayShootingStrategy(),
                       strategy])
        ensA = self.network.sampling_transitions[0].ensembles[0]
        ensA_sig = ((ensA,),(ensA,))
        strategy.group_weights['repex'] = 3.0
        strategy.mover_weights[('shooting',ensA_sig)] = 2.0

        root = scheme.move_decision_tree()
        (old_group_weights, old_mover_weights) = strategy.get_weights(
            scheme=scheme,
            sorted_movers=scheme.movers,
            sort_weights_override=strategy.group_weights
        )

        old_choice_prob = scheme.choice_probability
        old_movers = scheme.movers
        strategy.group_weights = {}
        strategy.mover_weights = {}

        (group_weights, mover_weights) = strategy.get_weights(
            scheme=scheme,
            sorted_movers=scheme.movers,
            sort_weights_override=strategy.group_weights,
            mover_weights_override=strategy.mover_weights
        )
        scheme.choice_probability = {}
        new_choice_prob = strategy.choice_probability(
            scheme, group_weights, mover_weights
        )

        for group in group_weights:
            assert_almost_equal(group_weights[group],
                                old_group_weights[group])

        for group in scheme.movers:
            for mover in scheme.movers[group]:
                old_mover = [
                    old for old in old_movers[group] 
                    if mover.ensemble_signature==old.ensemble_signature
                ][0]
                assert_almost_equal(new_choice_prob[mover], 
                                    old_choice_prob[old_mover])

        #print new_choice_prob
        for (old, new) in zip(list(old_mover_weights.keys()), list(mover_weights.keys())):
            try:
                assert_equal(old_mover_weights[old], mover_weights[new])
            except AssertionError:
                print(old_mover_weights[old])
                print(mover_weights[new])
                raise
        assert_equal(old_mover_weights, mover_weights)

    def test_get_weights_mover_weights_set_no_shooting(self):
        # follows test_get_weights_mover_weights_set, replacing shooting
        # with path reversal
        strategy = OrganizeByMoveGroupStrategy()
        scheme = MoveScheme(self.network)
        scheme.append([NearestNeighborRepExStrategy(), 
                       PathReversalStrategy(),
                       strategy])
        strategy.group_weights['repex'] = 3.0
        strategy.group_weights['pathreversal'] = 1.0
        root = scheme.move_decision_tree()
        (group_weights, mover_weights) = strategy.get_weights(
            scheme=scheme,
            sorted_movers=scheme.movers,
            sort_weights_override=strategy.group_weights
        )
        assert_equal(group_weights, {'pathreversal' : 1.0, 'repex' : 3.0})
        ensA = self.network.sampling_transitions[0].ensembles[0]
        ensA_sig = ((ensA,),(ensA,)) 
        strategy.group_weights = {}
        strategy.mover_weights[('pathreversal',ensA_sig)] = 2.0

        (group_weights, mover_weights) = strategy.get_weights(
            scheme=scheme,
            sorted_movers=scheme.movers,
            sort_weights_override=strategy.group_weights,
            mover_weights_override=strategy.mover_weights
        )
        for sig in [s for s in mover_weights if s[0]=='pathreversal']:
            weight_sigA = mover_weights[('pathreversal',ensA_sig)]
            ratio = old_div(mover_weights[sig], weight_sigA)
            if sig == ('pathreversal',ensA_sig):
                assert_equal(ratio, 1.0)
            else:
                assert_equal(ratio, 0.5)

        assert_almost_equal(group_weights['pathreversal'], 1.0)
        assert_almost_equal(group_weights['repex'], 3.0)


class TestOrganizeByEnsembleStrategy(MoveStrategyTestSetup):
    StrategyClass = OrganizeByEnsembleStrategy
    def setup(self):
        super(TestOrganizeByEnsembleStrategy, self).setup()
        scheme = MoveScheme(self.network)
        ens0 = self.network.sampling_transitions[0].ensembles[0]
        ens1 = self.network.sampling_transitions[0].ensembles[1]
        ens2 = self.network.sampling_transitions[0].ensembles[2]
        scheme.movers['shooting'] = [
            paths.OneWayShootingMover(
                selector=paths.UniformSelector(),
                ensemble=ens
            )
            for ens in [ens0, ens1, ens2]
        ]
        scheme.movers['repex'] = [
            paths.ReplicaExchangeMover(ens0, ens1),
            paths.ReplicaExchangeMover(ens1, ens2)
        ]
        scheme.movers['pathreversal'] = [
            paths.PathReversalMover(ensemble=ens) 
            for ens in [ens0, ens1, ens2]
        ]
        scheme.movers['minus'] = [paths.MinusMover(
            minus_ensemble=self.network.minus_ensembles[0],
            innermost_ensembles=[ens0]
        )]
        self.scheme = scheme

    def test_choice_probability(self):
        scheme = self.scheme
        ens0 = self.network.sampling_transitions[0].ensembles[0]
        ens1 = self.network.sampling_transitions[0].ensembles[1]
        ens2 = self.network.sampling_transitions[0].ensembles[2]
        minus = self.network.minus_ensembles[0]
        ens0_sig = ((ens0,),(ens0,))
        ensemble_weights = {ens0 : 2.0, ens1 : 1.0, ens2 : 1.0, minus : 0.5 }
        mover_weights = {}
        for groupname in scheme.movers:
            for m in scheme.movers[groupname]:
                sig = m.ensemble_signature
                for e in sig[0]:
                    mover_weights[(groupname, sig, e)] = 1.0
        mover_weights[('shooting', ens0_sig, ens0)] = 2.0
        strategy = self.StrategyClass()

        choice_prob = strategy.choice_probability(scheme, ensemble_weights,
                                                  mover_weights)
        # choice_prob is {mover : prob}; found is {mover_sig : prob}
        found = {
            (
                [g for g in scheme.movers if m in scheme.movers[g]][0],
                m.ensemble_signature
            ) : choice_prob[m] for m in choice_prob
        }

        # norm for ensemble selection: 4.5
        # ens0: 2/4.5 = 4/9
        # ens1, ens2: 2/9
        # minus: 1/9

        # ens0 moves: shooting (2.0), repex01 (1.0), minus (1.0), rev (1.0)
        # ens1 moves: shooting (1.0), repex01 (1.0), repex12 (1.0), rev (1.0) 
        # ens2 moves: shooting (1.0), repex12 (1.0), rev (1.0) 
        # minus moves: minus (1.0)

        # ens0 norm: 5.0
        # ens1 norm: 4.0
        # ens2 norm: 3.0
        # ens3 norm: 1.0

        # Probabilities
        #   shooting0 = 4/9 * 2.0/5.0 = 8.0/45.0
        #   shooting1 = 2/9 * 1.0/4.0 = 1.0/18.0
        #   shooting2 = 2/9 * 1.0/3.0 = 2.0/27.0
        #   repex01 = 4/9 * 1.0/5.0 + 2/9 * 1.0/4.0 = 13.0/90.0 
        #   repex12 = 2/9 * 1.0/4.0 + 2/9 * 1.0/3.0 = 7.0/54.0
        #   minus = 1/9*1.0 + 4/9*1.0/5.0 = 1.8/9.0
        #   rev0 = 4/9 * 1.0/5.0 = 4.0/45.0
        #   rev1 = 2/9 * 1.0/4.0 = 1.0/18.0
        #   rev2 = 2/9 * 1.0/3.0 = 2.0/27.0
        found_sigs = set([s[1] for s in list(found.keys())])

        sig0 = reorder_ensemble_signature(((ens0,),(ens0,)), found_sigs)
        sig1 = reorder_ensemble_signature(((ens1,),(ens1,)), found_sigs)
        sig2 = reorder_ensemble_signature(((ens2,),(ens2,)), found_sigs)
        sig01 = reorder_ensemble_signature(((ens0,ens1),(ens0,ens1)),
                                           found_sigs)
        sig12 = reorder_ensemble_signature(((ens1,ens2),(ens1,ens2)),
                                           found_sigs)
        sig_minus = reorder_ensemble_signature(((minus,ens0),(minus,ens0)),
                                               found_sigs)

        expected = {
            ('shooting', sig0) : old_div(8.0,45.0),
            ('shooting', sig1) : old_div(1.0,18.0),
            ('shooting', sig2) : old_div(2.0,27.0),
            ('repex', sig01) : old_div(13.0,90.0),
            ('repex', sig12) : old_div(7.0,54.0),
            ('minus', sig_minus) : old_div(1.8,9.0),
            ('pathreversal', sig0) : old_div(4.0,45.0),
            ('pathreversal', sig1) : old_div(1.0,18.0),
            ('pathreversal', sig2) : old_div(2.0,27.0)
        }
        assert_equal(set(expected.keys()), set(found.keys()))
        for k in list(expected.keys()):
            assert_almost_equal(expected[k], found[k])


    def test_chooser_mover_weights(self):
        scheme = self.scheme
        ens0 = self.network.sampling_transitions[0].ensembles[0]
        ens1 = self.network.sampling_transitions[0].ensembles[1]
        ens2 = self.network.sampling_transitions[0].ensembles[2]
        minus = self.network.minus_ensembles[0]
        strategy = self.StrategyClass()

        (ensemble_weights, mover_weights) = strategy.default_weights(scheme)

        for ens in [ens0, ens1, ens2, minus]:
            chooser_mweights = strategy.chooser_mover_weights(scheme, ens, 
                                                              mover_weights)
            for mover in chooser_mweights:
                assert_in(ens, mover.ensemble_signature[0])

            if ens in [ens0, ens1]:
                assert_equal(len(chooser_mweights), 4)
            elif ens is minus:
                assert_equal(len(chooser_mweights), 1)
            elif ens is ens2:
                assert_equal(len(chooser_mweights), 3)
        # that test feels a little minimal, but I guess it does the job

    def test_default_weights(self):
        scheme = self.scheme
        ens0 = self.network.sampling_transitions[0].ensembles[0]
        ens1 = self.network.sampling_transitions[0].ensembles[1]
        ens2 = self.network.sampling_transitions[0].ensembles[2]
        minus = self.network.minus_ensembles[0]
        strategy = self.StrategyClass()

        (ensemble_weights, mover_weights)= strategy.default_weights(scheme)

        assert_equal(ensemble_weights, 
                     {e : 1.0 for e in [ens0, ens1, ens2, minus]})

        # get the correct order on repex and minus signatures for testing
        double_sigs = [m.ensemble_signature 
                       for m in sum([scheme.movers[g] 
                                     for g in ['repex', 'minus']], [])
                      ]
        repex01_sig = reorder_ensemble_signature(((ens0,ens1),(ens0,ens1)),
                                                 double_sigs)
        repex12_sig = reorder_ensemble_signature(((ens1,ens2),(ens1,ens2)),
                                                 double_sigs)
        minus_sig = reorder_ensemble_signature(((minus,ens0),(minus,ens0)),
                                               double_sigs)

        mover_ens_sigs = [
            ('shooting', ((ens0,),(ens0,)), ens0),
            ('shooting', ((ens1,),(ens1,)), ens1),
            ('shooting', ((ens2,),(ens2,)), ens2),
            ('pathreversal', ((ens0,),(ens0,)), ens0),
            ('pathreversal', ((ens1,),(ens1,)), ens1),
            ('pathreversal', ((ens2,),(ens2,)), ens2),
            ('repex', repex01_sig, ens0),
            ('repex', repex01_sig, ens1),
            ('repex', repex12_sig, ens1),
            ('repex', repex12_sig, ens2),
            ('minus', minus_sig, minus),
            ('minus', minus_sig, ens0),
        ]

        assert_equal(len(mover_ens_sigs), len(list(mover_weights.keys())))
        assert_equal(set(mover_weights.keys()), set(mover_ens_sigs))
        
        expected = {s : 1.0 for s in mover_ens_sigs}
        assert_equal(expected, mover_weights)

    def test_weights_from_choice_probability(self):
        scheme = self.scheme
        ens0 = self.network.sampling_transitions[0].ensembles[0]
        ens1 = self.network.sampling_transitions[0].ensembles[1]
        ens2 = self.network.sampling_transitions[0].ensembles[2]
        minus = self.network.minus_ensembles[0]

        sig0 = ((ens0,),(ens0,))
        sig1 = ((ens1,),(ens1,))
        sig2 = ((ens2,),(ens2,))
        sig01 = ((ens0,ens1),(ens0,ens1))
        sig12 = ((ens1,ens2),(ens1,ens2))
        sig_minus = ((minus,ens0),(minus,ens0))

        choice_probability = {
            find_mover(scheme, 'shooting', sig0) : old_div(8.0,45.0),
            find_mover(scheme, 'shooting', sig1) : old_div(1.0,18.0),
            find_mover(scheme, 'shooting', sig2) : old_div(2.0,27.0),
            find_mover(scheme, 'repex', sig01) : old_div(13.0,90.0),
            find_mover(scheme, 'repex', sig12) : old_div(7.0,54.0),
            find_mover(scheme, 'minus', sig_minus) : old_div(1.8,9.0),
            find_mover(scheme, 'pathreversal', sig0) : old_div(4.0,45.0),
            find_mover(scheme, 'pathreversal', sig1) : old_div(1.0,18.0),
            find_mover(scheme, 'pathreversal', sig2) : old_div(2.0,27.0)
        }

        strategy = self.StrategyClass()
        (ens_w, mover_w) = strategy.weights_from_choice_probability(
            scheme, choice_probability
        )

        new_choice_prob = strategy.choice_probability(scheme, ens_w, mover_w)
        assert_equal_=(choice_probability, new_choice_prob)


    def test_make_movers(self):
        scheme = self.scheme
        strategy = self.StrategyClass()
        root = strategy.make_movers(scheme)

        choosers = root.movers
        chooser_weights = root.weights
        assert_equal(list(chooser_weights), [1.0]*4)
        for (mover, w) in zip(choosers, chooser_weights):
            n_moves = len(mover.movers)
            assert_equal(mover.weights, [1.0]*n_moves)


    def test_make_mover_rebuild_choice_probability(self):
        scheme = self.scheme

        # Organize by move group, switch to ensemble, switch back
        scheme.append(OrganizeByMoveGroupStrategy())
        root_1a = scheme.move_decision_tree()
        choice_prob_1a = scheme.choice_probability
        scheme.append(self.StrategyClass(), force=True)
        root_1b = scheme.move_decision_tree(rebuild=True)
        choice_prob_1b = scheme.choice_probability
        assert(choice_prob_1a is not choice_prob_1b)
        for m in choice_prob_1a:
            assert_almost_equal(choice_prob_1a[m], choice_prob_1b[m])
        scheme.append(OrganizeByMoveGroupStrategy(), force=True)
        root_1c = scheme.move_decision_tree(rebuild=True)
        choice_prob_1c = scheme.choice_probability
        assert(choice_prob_1a is not choice_prob_1c)
        for m in choice_prob_1a:
            assert_almost_equal(choice_prob_1a[m], choice_prob_1c[m])

        # Organize by ensemble, switch to move group, switch back
        scheme.strategies.clear()
        scheme.append(self.StrategyClass(), force=True)
        root_2a = scheme.move_decision_tree(rebuild=True)
        choice_prob_2a = scheme.choice_probability
        scheme.append(OrganizeByMoveGroupStrategy(), force=True)
        root_2b = scheme.move_decision_tree(rebuild=True)
        choice_prob_2b = scheme.choice_probability
        assert(choice_prob_2a is not choice_prob_2b)
        for m in choice_prob_2a:
            assert_almost_equal(choice_prob_2a[m], choice_prob_2b[m])
        scheme.append(self.StrategyClass(), force=True)
        root_2c = scheme.move_decision_tree(rebuild=True)
        choice_prob_2c = scheme.choice_probability
        assert(choice_prob_1a is not choice_prob_1c)
        for m in choice_prob_2a:
            assert_almost_equal(choice_prob_2a[m], choice_prob_2c[m])

        # Org strategy 1 and org strategy 2 are different
        for m in choice_prob_1a:
            assert(abs(choice_prob_1a[m] - choice_prob_2a[m]) > 0.001)

class TestPoorSingleReplicaStrategy(TestOrganizeByEnsembleStrategy):
    StrategyClass = PoorSingleReplicaStrategy

    def test_chooser_mover_weights(self):
        scheme = self.scheme
        ens0 = self.network.sampling_transitions[0].ensembles[0]
        ens1 = self.network.sampling_transitions[0].ensembles[1]
        ens2 = self.network.sampling_transitions[0].ensembles[2]
        minus = self.network.minus_ensembles[0]
        strategy = self.StrategyClass()

        (ensemble_weights, mover_weights) = strategy.default_weights(scheme)
        (ens_w, mov_w) = strategy.get_weights(scheme, scheme.movers,
                                              strategy.ensemble_weights,
                                              strategy.mover_weights)
        scheme.choice_probability = strategy.choice_probability(scheme, 
                                                                ens_w, mov_w)
        for ens in [ens0, ens1, ens2, minus]:
            chooser_mweights = strategy.chooser_mover_weights(scheme, ens, 
                                                              mover_weights)
            if ens in [ens0, ens1]:
                assert_equal(len(chooser_mweights), 5)
            elif ens is minus:
                assert_equal(len(chooser_mweights), 2)
            elif ens is ens2:
                assert_equal(len(chooser_mweights), 4)

            real_movers = [m for m in list(chooser_mweights.keys()) 
                           if m != strategy.null_mover]
            for m in real_movers:
                assert_equal(scheme.choice_probability[m], chooser_mweights[m])
            assert_almost_equal(sum(chooser_mweights.values()), 1.0)

    def test_make_movers(self):
        scheme = self.scheme
        strategy = self.StrategyClass()
        root = strategy.make_movers(scheme)

        assert_equal(len(root.movers), 4)

        for mover in list(scheme.choice_probability.keys()):
            assert_almost_equal(
                scheme.choice_probability[mover] * 0.25,
                scheme.real_choice_probability[mover]
            )



