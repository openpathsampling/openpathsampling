from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        assert_almost_equal, raises, assert_in)
from nose.plugins.skip import Skip, SkipTest
from test_helpers import true_func, assert_equal_array_array, make_1d_traj

import openpathsampling as paths
from openpathsampling.analysis.move_scheme import MoveScheme, DefaultScheme
from openpathsampling.analysis.move_strategy import *
from openpathsampling import VolumeFactory as vf

import collections

import logging
logging.getLogger('openpathsampling.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.ensemble').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)


class testStrategyLevels(object):
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
        cvA = paths.CV_Function(name="xA", fcn=lambda s : s.xyz[0][0])
        cvB = paths.CV_Function(name="xB", fcn=lambda s : -s.xyz[0][0])
        self.stateA = paths.CVRangeVolume(cvA, float("-inf"), -0.5)
        self.stateB = paths.CVRangeVolume(cvB, float("-inf"), -0.5)
        interfacesA = vf.CVRangeVolumeSet(cvA, float("-inf"), 
                                          [-0.5, -0.3, -0.1, 0.0])
        interfacesB = vf.CVRangeVolumeSet(cvB, float("-inf"), 
                                          [-0.5, -0.3, -0.1, 0.0])
        self.network = paths.MSTISNetwork([
            (self.stateA, interfacesA, "A", cvA),
            (self.stateB, interfacesB, "B", cvB)
        ])


class testMoveStrategy(MoveStrategyTestSetup):
    def test_levels(self):
        strategy = MoveStrategy(ensembles=None, group="test", replace=True)
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
        self.strategy = MoveStrategy(ensembles=None, group="test", replace=True)
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

class testOneWayShootingStrategy(MoveStrategyTestSetup):
    def test_make_movers(self):
        strategy = OneWayShootingStrategy()
        scheme = MoveScheme(self.network)
        movers = strategy.make_movers(scheme)
        assert_equal(len(movers), 6)
        for mover in movers:
            assert_equal(type(mover), paths.OneWayShootingMover)
            assert_equal(type(mover.selector), paths.UniformSelector)

class testNearestNeighborRepExStrategy(MoveStrategyTestSetup):
    def test_make_movers(self):
        strategy = NearestNeighborRepExStrategy()
        scheme = MoveScheme(self.network)
        movers = strategy.make_movers(scheme)
        assert_equal(len(movers), 4)
        ens0 = self.network.sampling_transitions[0].ensembles
        ens1 = self.network.sampling_transitions[1].ensembles
        assert_equal_array_array(movers[0].ensemble_signature,
                     ((ens0[0], ens0[1]), (ens0[0], ens0[1])))
        assert_equal_array_array(movers[1].ensemble_signature,
                     ((ens0[1], ens0[2]), (ens0[1], ens0[2])))
        assert_equal_array_array(movers[2].ensemble_signature,
                     ((ens1[0], ens1[1]), (ens1[0], ens1[1])))
        assert_equal_array_array(movers[3].ensemble_signature,
                     ((ens1[1], ens1[2]), (ens1[1], ens1[2])))

class testAllSetRepExStrategy(MoveStrategyTestSetup):
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

class testSelectedPairsRepExStrategy(MoveStrategyTestSetup):
    def test_make_movers(self):
        ens00 = self.network.sampling_transitions[0].ensembles[0]
        ens02 = self.network.sampling_transitions[0].ensembles[2]
        strategy = SelectedPairsRepExStrategy(ensembles=[ens00, ens02])
        scheme = MoveScheme(self.network)
        movers = strategy.make_movers(scheme)
        assert_equal(len(movers), 1)
        assert_equal(movers[0].ensemble_signature_set, 
                     (set([ens00, ens02]), (set([ens00, ens02]))))

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
                     (set([ens00, ens01]), set([ens00, ens01])))
        assert_equal(movers[1].ensemble_signature_set,
                     (set([ens00, ens02]), set([ens00, ens02])))
        assert_equal(movers[2].ensemble_signature_set,
                     (set([ens01, ens02]), set([ens01, ens02])))


class testPathReversalStrategy(MoveStrategyTestSetup):
    def test_make_movers(self):
        strategy = PathReversalStrategy()
        scheme = MoveScheme(self.network)
        movers = strategy.make_movers(scheme)
        assert_equal(len(movers), 6)
        for m in movers:
            assert_equal(type(m), paths.PathReversalMover)


class testMinusMoveStrategy(MoveStrategyTestSetup):
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
        extra_minus = paths.MinusInterfaceEnsemble(
            state_vol=self.network.sampling_transitions[0].stateA,
            innermost_vols=[innerA, innerB]
        )
        self.network.special_ensembles['minus'][extra_minus] = [innerA, innerB]
        scheme = MoveScheme(self.network)
        ensembles = strategy.get_ensembles(scheme, None)
        assert_equal(len(ensembles), 2)
        assert_equal(set([len(ensembles[0]), len(ensembles[1])]), set([1,2]))

    def test_get_ensembles_fixed_ensembles(self):
        strategy = MinusMoveStrategy()
        minusA = self.network.special_ensembles['minus'].keys()[0]
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
        ens_minusA = minuses.keys()[0]
        ens_innerA = [t.ensembles[0] for t in minuses[ens_minusA]]
        sig_A = set([ens_minusA] + ens_innerA)
        ens_minusB = minuses.keys()[1]
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


class testOrganizeByEnsembleStrategy(MoveStrategyTestSetup):
    def setup(self):
        super(testOrganizeByEnsembleStrategy, self).setup()
        scheme = MoveScheme(self.network)
        scheme.movers = {} # handles LEGACY stuff
        ens0 = self.network.sampling_transitions[0].ensembles[0]
        ens1 = self.network.sampling_transitions[0].ensembles[1]
        ens2 = self.network.sampling_transitions[0].ensembles[2]
        scheme.movers['shooting'] = [
            paths.OneWayShootingMover(
                selector=paths.UniformSelector(),
                ensembles=[ens]
            )
            for ens in [ens0, ens1, ens2]
        ]
        scheme.movers['repex'] = [
            paths.ReplicaExchangeMover(ensembles=[ens0, ens1]),
            paths.ReplicaExchangeMover(ensembles=[ens1, ens2])
        ]
        scheme.movers['pathreversal'] = [
            paths.PathReversalMover(ensembles=ens) 
            for ens in [ens0, ens1, ens2]
        ]
        scheme.movers['minus'] = [paths.MinusMover(
            minus_ensemble=self.network.minus_ensembles[0],
            innermost_ensembles=[ens0]
        )]
        self.scheme = scheme

    def test_make_movers(self):
        scheme = self.scheme
        strategy = OrganizeByEnsembleStrategy()
        root = strategy.make_movers(scheme)

        choosers = root.movers
        chooser_weights = root.weights
        assert_equal(collections.Counter(chooser_weights), 
                     collections.Counter([0.5, 2.5, 3.0, 3.0]))
        for (mover, w) in zip(choosers, chooser_weights):
            weights_count = collections.Counter(mover.weights)

            if len(mover.weights) == 1: # minus
                assert_equal(weights_count, collections.Counter({ 0.5 : 1}))
                assert_equal(w, 0.5)
            elif len(mover.weights) == 3: # outermost
                assert_equal(weights_count, 
                             collections.Counter({ 0.5 : 1, 1.0 : 2}))
                assert_equal(w, 2.5)
            elif len(mover.weights) == 4:
                assert_equal(weights_count, 
                             collections.Counter({ 0.5 : 2, 1.0 : 2}))
                assert_equal(w, 3.0)

        for mover in scheme.choice_probability:
            assert_almost_equal(scheme.choice_probability[mover], 1.0/9.0)

    def test_make_movers_light_minus_ensemble(self):
        scheme = self.scheme
        strategy = OrganizeByEnsembleStrategy()
        minus =  scheme.network.special_ensembles['minus'].keys()[0]
        strategy.ensemble_weights[minus] = 0.5
        root = strategy.make_movers(scheme)

        choosers = root.movers
        chooser_weights = root.weights
        assert_equal(collections.Counter(chooser_weights),
                     collections.Counter([0.25, 2.5, 3.0, 3.0]))
        for mover in scheme.choice_probability:
            if minus in mover.input_ensembles:
                assert_equal(scheme.choice_probability[mover], 0.75/8.75)
            else:
                assert_equal(scheme.choice_probability[mover], 1.0/8.75)


    def test_make_mover_light_minus_and_minus_partner(self):
        scheme = self.scheme
        strategy = OrganizeByEnsembleStrategy()
        minus_mover = scheme.movers['minus'][0]
        minus_ensemble = minus_mover.minus_ensemble
        innermost_ensemble = minus_mover.innermost_ensemble
        minus_key = strategy._mover_key(scheme.movers['minus'][0], scheme)
        # this has no effect: minus mover is only move in minus ens
        strategy.mover_weights[minus_ensemble] = { minus_key : 0.5 }
        strategy.mover_weights[innermost_ensemble] = { minus_key : 0.5 }
        # note that if you really want to reduce the probability of the
        # specific move, it might be best to do that within
        # scheme.choice_probability

        root = strategy.make_movers(scheme)

        choosers = root.movers
        chooser_weights = root.weights
        for (chooser, w) in zip(choosers, chooser_weights):
            assert_almost_equal(sum(chooser.weights), w)

        assert_equal(collections.Counter(chooser_weights),
                     collections.Counter([0.5, 2.5, 2.75, 3.0]))

        for mover in scheme.choice_probability:
            if strategy._mover_key(mover, scheme) == minus_key:
                assert_almost_equal(scheme.choice_probability[mover], 0.75/8.75)
            else:
                assert_almost_equal(scheme.choice_probability[mover], 1.0/8.75)

    def test_make_mover_rebuild_choice_probability(self):
        raise SkipTest

    def test_make_movers_preserves_default_choice_prob(self):
        scheme1 = DefaultScheme(self.network)
        scheme1.movers = {} # handles LEGACY stuff
        scheme1.append(OrganizeByEnsembleStrategy())
        root1 = scheme1.move_decision_tree()

        scheme2 = DefaultScheme(self.network)
        root2 = scheme2.move_decision_tree()

        print scheme1.choice_probability

        raise SkipTest
                

class testDefaultStrategy(MoveStrategyTestSetup):
    def test_make_movers(self):
        scheme = MoveScheme(self.network)
        scheme.movers = {} # handles LEGACY stuff
        ens0 = self.network.sampling_transitions[0].ensembles[0]
        ens1 = self.network.sampling_transitions[0].ensembles[1]
        ens2 = self.network.sampling_transitions[0].ensembles[2]
        scheme.movers['shooting'] = [
            paths.OneWayShootingMover(
                selector=paths.UniformSelector(),
                ensembles=[ens]
            )
            for ens in [ens0, ens1, ens2]
        ]
        scheme.movers['repex'] = [
            paths.ReplicaExchangeMover(ensembles=[ens0, ens1]),
            paths.ReplicaExchangeMover(ensembles=[ens1, ens2])
        ]
        scheme.movers['pathreversal'] = [
            paths.PathReversalMover(ensembles=ens) 
            for ens in [ens0, ens1, ens2]
        ]
        scheme.movers['minus'] = [paths.MinusMover(
            minus_ensemble=self.network.minus_ensembles[0],
            innermost_ensembles=[ens0]
        )]

        strategy = DefaultStrategy()
        root = strategy.make_movers(scheme)
        
        assert_equal(len(root.movers), 4)
        names = ['ShootingChooser', 'RepexChooser', 'PathreversalChooser', 
                 'MinusChooser']
        name_dict = {root.movers[i].name : i for i in range(len(root.movers))}
        for name in names:
            assert_in(name, name_dict.keys())

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
        scheme.movers = {} # handles LEGACY stuff
        ens0 = self.network.sampling_transitions[0].ensembles[0]
        ens1 = self.network.sampling_transitions[0].ensembles[1]
        ens2 = self.network.sampling_transitions[0].ensembles[2]
        scheme.movers['blahblah']  = [
            paths.ReplicaExchangeMover(ensembles=[ens0, ens1]),
            paths.ReplicaExchangeMover(ensembles=[ens1, ens2])
        ]

        strategy = DefaultStrategy()
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
        scheme.movers = {} # handles LEGACY stuff
        ens0 = self.network.sampling_transitions[0].ensembles[0]
        ens1 = self.network.sampling_transitions[0].ensembles[1]
        ens2 = self.network.sampling_transitions[0].ensembles[2]
        scheme.movers['blahblahblah']  = [
            paths.ReplicaExchangeMover(ensembles=[ens0, ens1]),
            paths.ReplicaExchangeMover(ensembles=[ens1, ens2])
        ]

        strategy = DefaultStrategy()
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
        strategy = DefaultStrategy()

        scheme = MoveScheme(self.network)
        scheme.movers = {} # handles LEGACY stuff
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
            sorted_movers=scheme.movers,
            preset_sortkey_weights=strategy.group_weights
        )
        assert_equal(group_weights, {'shooting' : 1.0, 'repex' : 0.5})
        for group in mover_weights:
            for sig in mover_weights[group]:
                assert_equal(mover_weights[group][sig], 1.0)
                assert_in(sig, all_movers_sigs)
            assert_equal(len(mover_weights[group]), len(scheme.movers[group]))

        # check that we can reuse these in a different scheme
        scheme2 = MoveScheme(self.network)
        scheme2.movers = {} # handles LEGACY stuff
        scheme2.append(OneWayShootingStrategy())
        root = scheme2.move_decision_tree()
        assert_equal(len(scheme2.movers), 1)

        (group_weights, mover_weights) = strategy.get_weights(
            scheme=scheme2,
            sorted_movers=scheme2.movers,
            preset_sortkey_weights=strategy.group_weights
        )
        assert_equal(group_weights, {'shooting' : 1.0})
        for group in mover_weights:
            for sig in mover_weights[group]:
                assert_equal(mover_weights[group][sig], 1.0)
                assert_in(sig, [m.ensemble_signature 
                                for m in scheme2.movers[group]])
            assert_equal(len(mover_weights[group]), len(scheme2.movers[group]))

    def test_get_weights_both_internal_weights_set(self):
        strategy = DefaultStrategy()
        ensA = self.network.sampling_transitions[0].ensembles[0]
        ensA_sig = ((ensA,),(ensA,))

        scheme = MoveScheme(self.network)
        scheme.movers = {} # handles LEGACY stuff
        scheme.append([NearestNeighborRepExStrategy(), 
                       OneWayShootingStrategy(), 
                       strategy])
        root = scheme.move_decision_tree()
        old_choice_probability = scheme.choice_probability
        old_shooter_ensA = [m for m in scheme.movers['shooting'] 
                            if m.ensemble_signature == ensA_sig][0]

        strategy.group_weights['repex'] = 3.0
        strategy.mover_weights['shooting'] = {} #sadly, can't safely avoid 
        strategy.mover_weights['shooting'][ensA_sig]= 2.0

        (group_weights, mover_weights) = strategy.get_weights(
            scheme=scheme,
            sorted_movers=scheme.movers,
            preset_sortkey_weights=strategy.group_weights
        )
        root = scheme.move_decision_tree(rebuild=True)

        assert_equal(group_weights, {'shooting' : 1.0, 'repex' : 3.0})
        expected_mover_weights = {}
        for group in scheme.movers:
            expected_mover_weights[group] = {}
            for mover in scheme.movers[group]:
                sig = mover.ensemble_signature
                if group == 'shooting' and sig == ensA_sig:
                    expected_mover_weights[group][sig] = 1.0
                elif group == 'shooting':
                    expected_mover_weights[group][sig] = 0.5
                else:
                    expected_mover_weights[group][sig] = 1.0

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
        assert_equal(new_shoot_chooser.weights[shooter_ensA_idx], 1.0)

        assert_not_equal(new_choice_probability, old_choice_probability)

    def test_get_weights_group_weights_set(self):
        strategy = DefaultStrategy()
        scheme = MoveScheme(self.network)
        scheme.movers = {} # handles LEGACY stuff
        scheme.append([NearestNeighborRepExStrategy(), 
                       OneWayShootingStrategy(), 
                       strategy])
        root = scheme.move_decision_tree()

        nshoot = len(scheme.movers['shooting'])
        nrepex = len(scheme.movers['repex'])

        for mover in scheme.movers['shooting']:
            assert_almost_equal(scheme.choice_probability[mover], 1.0/9.0)
        for mover in scheme.movers['repex']:
            assert_almost_equal(scheme.choice_probability[mover], 1.0/12.0)


        strategy.group_weights['shooting'] = 2.0
        #root = scheme.move_decision_tree(rebuild=True)

        assert_equal(strategy.group_weights, {'shooting' : 2.0})
        assert_equal(strategy.mover_weights, {})

        (group_weights, mover_weights) = strategy.get_weights(
            scheme=scheme,
            sorted_movers=scheme.movers,
            preset_sortkey_weights=strategy.group_weights
        )
        assert_equal(group_weights, {'shooting' : 2.0, 'repex' : 0.5})
        for group in mover_weights:
            for mover in mover_weights[group]:
                assert_equal(mover_weights[group][mover],
                             mover_weights[group].values()[0])

        choice_prob = strategy.choice_probability(
            scheme, scheme.movers, group_weights, mover_weights
        )
        for mover in scheme.movers['shooting']:
            assert_almost_equal(choice_prob[mover], 2.0/15.0)
        for mover in scheme.movers['repex']:
            assert_almost_equal(choice_prob[mover], 1.0/20.0)
        

    def test_get_weights_mover_weights_set(self):
        strategy = DefaultStrategy()
        scheme = MoveScheme(self.network)
        scheme.movers = {} # handles LEGACY stuff
        scheme.append([NearestNeighborRepExStrategy(), 
                       OneWayShootingStrategy(),
                       strategy])
        # build it so there's an old version
        strategy.group_weights['repex'] = 3.0
        root = scheme.move_decision_tree()
        (group_weights, mover_weights) = strategy.get_weights(
            scheme=scheme,
            sorted_movers=scheme.movers,
            preset_sortkey_weights=strategy.group_weights
        )
        assert_equal(group_weights, {'shooting' : 1.0, 'repex' : 3.0})

        ensA = self.network.sampling_transitions[0].ensembles[0]
        ensA_sig = ((ensA,),(ensA,)) 
        strategy.group_weights = {}
        strategy.mover_weights['shooting'] = {} #sadly, can't safely avoid 
        strategy.mover_weights['shooting'][ensA_sig] = 2.0

        (group_weights, mover_weights) = strategy.get_weights(
            scheme=scheme,
            sorted_movers=scheme.movers,
            preset_sortkey_weights=strategy.group_weights
        )
        #root = scheme.move_decision_tree(rebuild=True) 

        # we always normalize to so that shooting for ensA is 1.0
        for sig in mover_weights['shooting']:
            if sig == ensA_sig:
                assert_equal(mover_weights['shooting'][sig], 1.0)
            else:
                assert_equal(mover_weights['shooting'][sig], 0.5)

        assert_almost_equal(group_weights['shooting'], 1.0)
        assert_almost_equal(group_weights['repex'], 3.0)

    def test_get_weights_internal_unset_choice_prob_set(self):
        strategy = DefaultStrategy()
        scheme = MoveScheme(self.network)
        scheme.movers = {} # handles LEGACY stuff
        scheme.append([NearestNeighborRepExStrategy(), 
                       OneWayShootingStrategy(),
                       strategy])
        ensA = self.network.sampling_transitions[0].ensembles[0]
        ensA_sig = ((ensA,),(ensA,))
        strategy.group_weights['repex'] = 3.0
        strategy.mover_weights['shooting'] = {} #sadly, can't safely avoid 
        strategy.mover_weights['shooting'][ensA_sig] = 2.0

        root = scheme.move_decision_tree()
        (old_group_weights, old_mover_weights) = strategy.get_weights(
            scheme=scheme,
            sorted_movers=scheme.movers,
            preset_sortkey_weights=strategy.group_weights
        )

        old_choice_prob = scheme.choice_probability
        old_movers = scheme.movers
        strategy.group_weights = {}
        strategy.mover_weights = {}

        (group_weights, mover_weights) = strategy.get_weights(
            scheme=scheme,
            sorted_movers=scheme.movers,
            preset_sortkey_weights=strategy.group_weights
        )
        scheme.choice_probability = {}
        new_choice_prob = strategy.choice_probability(
            scheme, scheme.movers, group_weights, mover_weights
        )

        for group in scheme.movers:
            for mover in scheme.movers[group]:
                old_mover = [
                    old for old in old_movers[group] 
                    if mover.ensemble_signature==old.ensemble_signature
                ][0]
                assert_almost_equal(new_choice_prob[mover], 
                                    old_choice_prob[old_mover])

        for group in group_weights:
            assert_almost_equal(group_weights[group],
                                old_group_weights[group])

        #print new_choice_prob
        for (old, new) in zip(old_mover_weights.keys(), mover_weights.keys()):
            print old, new
            try:
                assert_equal(old_mover_weights[old], mover_weights[new])
        
            except AssertionError:
                print old_mover_weights[old]
                print mover_weights[new]
                raise
        assert_equal(old_mover_weights, mover_weights)

    def test_get_weights_mover_weights_set_no_shooting(self):
        # follows test_get_weights_mover_weights_set, replacing shooting
        # with path reversal
        strategy = DefaultStrategy()
        scheme = MoveScheme(self.network)
        scheme.movers = {} # handles LEGACY stuff
        scheme.append([NearestNeighborRepExStrategy(), 
                       PathReversalStrategy(),
                       strategy])
        strategy.group_weights['repex'] = 3.0
        strategy.group_weights['pathreversal'] = 1.0
        root = scheme.move_decision_tree()
        (group_weights, mover_weights) = strategy.get_weights(
            scheme=scheme,
            sorted_movers=scheme.movers,
            preset_sortkey_weights=strategy.group_weights
        )
        assert_equal(group_weights, {'pathreversal' : 1.0, 'repex' : 3.0})
        ensA = self.network.sampling_transitions[0].ensembles[0]
        ensA_sig = ((ensA,),(ensA,)) 
        strategy.group_weights = {}
        strategy.mover_weights['pathreversal'] = {} #sadly, can't safely avoid 
        strategy.mover_weights['pathreversal'][ensA_sig] = 2.0

        (group_weights, mover_weights) = strategy.get_weights(
            scheme=scheme,
            sorted_movers=scheme.movers,
            preset_sortkey_weights=strategy.group_weights
        )
        for sig in mover_weights['pathreversal']:
            weight_sigA = mover_weights['pathreversal'][ensA_sig]
            ratio = mover_weights['pathreversal'][sig] / weight_sigA
            if sig == ensA_sig:
                assert_equal(ratio, 1.0)
            else:
                assert_equal(ratio, 0.5)

        assert_almost_equal(group_weights['pathreversal'], 0.25)
        assert_almost_equal(group_weights['repex'], 0.75)

