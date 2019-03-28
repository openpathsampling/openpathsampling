from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from builtins import zip
from builtins import range
from builtins import object
from past.utils import old_div
from nose.tools import (assert_equal, assert_not_equal, assert_almost_equal,
                        assert_in, raises, assert_is, assert_is_not,
                        assert_true)
from nose.plugins.skip import Skip, SkipTest
from .test_helpers import (true_func, assert_equal_array_array,
                           make_1d_traj, assert_items_equal)

import copy

import openpathsampling as paths
from openpathsampling.high_level.move_scheme import *
from openpathsampling.high_level.move_strategy import (
    levels,
    MoveStrategy, OneWayShootingStrategy, NearestNeighborRepExStrategy,
    OrganizeByMoveGroupStrategy, AllSetRepExStrategy
)

import openpathsampling.high_level.move_strategy as strategies

import logging
logging.getLogger('openpathsampling.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.ensemble').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.netcdfplus').setLevel(logging.CRITICAL)

class TestMoveAcceptanceAnalysis(object):
    def setup(self):
        pass

    def test_initialization(self):
        pass

    def test_add_steps(self):
        pass

    def test_no_move_keys(self):
        pass

    def test_n_total_trials(self):
        pass

    def test_select_movers_groupname(self):
        pass

    def test_select_movers_mover(self):
        pass

    def test_summary_data_groupname(self):
        pass

    def test_summary_data_mover(self):
        pass

    def test_format_as_text(self):
        pass



class TestMoveScheme(object):
    def setup(self):
        paths.InterfaceSet._reset()
        cvA = paths.FunctionCV(name="xA", f=lambda s : s.xyz[0][0])
        cvB = paths.FunctionCV(name="xB", f=lambda s : -s.xyz[0][0])
        self.stateA = paths.CVDefinedVolume(cvA, float("-inf"), -0.5)
        self.stateB = paths.CVDefinedVolume(cvB, float("-inf"), -0.5)
        interfacesA = paths.VolumeInterfaceSet(cvA, float("-inf"), 
                                               [-0.5, -0.3, -0.1])
        interfacesB = paths.VolumeInterfaceSet(cvB, float("-inf"), 
                                               [-0.5, -0.3, -0.1])
        network = paths.MSTISNetwork(
            [(self.stateA, interfacesA),
             (self.stateB, interfacesB)],
            ms_outers=paths.MSOuterTISInterface.from_lambdas(
                {interfacesA: 0.0, interfacesB: 0.0}
            )
        )
        self.scheme = MoveScheme(network)

    def test_append_individuals_default_levels(self):
        shootstrat = OneWayShootingStrategy()
        repexstrat = NearestNeighborRepExStrategy()
        defaultstrat = OrganizeByMoveGroupStrategy()
        assert_equal(len(list(self.scheme.strategies.keys())), 0)
        self.scheme.append(shootstrat)
        self.scheme.append(repexstrat)
        self.scheme.append(defaultstrat)

        strats = self.scheme.strategies
        assert_equal(len(list(strats.keys())), 3)
        pairs = [(levels.MOVER, shootstrat), (levels.SIGNATURE, repexstrat), 
                 (levels.GLOBAL, defaultstrat)]
        for (k, v) in pairs:
            assert_in(v, strats[k])

    def test_append_groups_default_levels(self):
        shootstrat = OneWayShootingStrategy()
        repexstrat = NearestNeighborRepExStrategy()
        defaultstrat = OrganizeByMoveGroupStrategy()
        assert_equal(len(list(self.scheme.strategies.keys())), 0)
        self.scheme.append([shootstrat, repexstrat, defaultstrat])

        strats = self.scheme.strategies
        assert_equal(len(list(strats.keys())), 3)
        pairs = [(levels.MOVER, shootstrat), (levels.SIGNATURE, repexstrat), 
                 (levels.GLOBAL, defaultstrat)]
        for (k, v) in pairs:
            assert_in(v, strats[k])


    def test_append_individuals_custom_levels(self):
        shootstrat = OneWayShootingStrategy()
        repexstrat = NearestNeighborRepExStrategy()
        defaultstrat = OrganizeByMoveGroupStrategy()
        assert_equal(len(list(self.scheme.strategies.keys())), 0)
        self.scheme.append(shootstrat, 60)
        self.scheme.append(repexstrat, 60)
        self.scheme.append(defaultstrat, 60)

        strats = self.scheme.strategies
        assert_equal(len(list(strats.keys())), 1)
        assert_items_equal(strats[60], [shootstrat, repexstrat, defaultstrat])

    def test_append_groups_same_custom_level(self):
        shootstrat = OneWayShootingStrategy()
        repexstrat = NearestNeighborRepExStrategy()
        defaultstrat = OrganizeByMoveGroupStrategy()
        assert_equal(len(list(self.scheme.strategies.keys())), 0)
        self.scheme.append([shootstrat, repexstrat, defaultstrat], 60)

        strats = self.scheme.strategies
        assert_equal(len(list(strats.keys())), 1)
        assert_items_equal(strats[60], [shootstrat, repexstrat, defaultstrat])

    def test_append_group_different_custom_levels(self):
        shootstrat = OneWayShootingStrategy()
        repexstrat = NearestNeighborRepExStrategy()
        defaultstrat = OrganizeByMoveGroupStrategy()
        assert_equal(len(list(self.scheme.strategies.keys())), 0)
        self.scheme.append([shootstrat, repexstrat, defaultstrat], 
                           [45, 55, 65])

        strats = self.scheme.strategies
        assert_equal(len(list(strats.keys())), 3)
        for (k, v) in [(45, shootstrat), (55, repexstrat), (65, defaultstrat)]:
            assert_in(v, strats[k])

    def test_apply_strategy(self):
        if self.scheme.movers == {}:
            print("Full support of MoveStrategy implemented?")
            print("Time to remove legacy from tests.")
        else:
            self.scheme.movers = {}

        shoot_strat_1 = OneWayShootingStrategy(
            ensembles=self.scheme.network.sampling_transitions[0].ensembles,
            replace=False
        )
        shoot_strat_2 = OneWayShootingStrategy(
            ensembles=(
                [self.scheme.network.sampling_transitions[0].ensembles[-1]] +
                self.scheme.network.sampling_transitions[1].ensembles
            ),
            replace=False
        )
        shoot_strat_3 = OneWayShootingStrategy(replace=True)

        self.scheme.apply_strategy(shoot_strat_1)
        assert_items_equal(list(self.scheme.movers.keys()), ['shooting'])
        assert_equal(len(self.scheme.movers['shooting']), 3)

        self.scheme.apply_strategy(shoot_strat_2)
        assert_items_equal(list(self.scheme.movers.keys()), ['shooting'])
        assert_equal(len(self.scheme.movers['shooting']), 7)
        old_movers = copy.copy(self.scheme.movers['shooting'])

        self.scheme.apply_strategy(shoot_strat_3)
        assert_items_equal(list(self.scheme.movers.keys()), ['shooting'])
        assert_equal(len(self.scheme.movers['shooting']), 7)
        new_movers = self.scheme.movers['shooting']
        for (o, n) in zip(old_movers, new_movers):
            assert_equal(o is n, False)

        shoot_strat_3.replace_signatures = True
        self.scheme.apply_strategy(shoot_strat_3)
        assert_equal(len(self.scheme.movers['shooting']), 6)

        self.scheme.movers = {}
        shoot_strat_1.set_replace(True)
        self.scheme.apply_strategy(shoot_strat_1)
        assert_equal(len(self.scheme.movers['shooting']), 3)
        old_movers = copy.copy(self.scheme.movers['shooting'])

        shoot_strat_3.replace_signatures = False
        self.scheme.apply_strategy(shoot_strat_3)
        assert_equal(len(self.scheme.movers['shooting']), 6)
        new_movers = self.scheme.movers['shooting']
        for (o, n) in zip(old_movers, new_movers):
            assert_equal(o is n, False)

    def test_move_decision_tree(self):
        shoot = OneWayShootingStrategy()
        repex = NearestNeighborRepExStrategy()
        default = OrganizeByMoveGroupStrategy()
        self.scheme.append([default, shoot, repex])

        assert_equal(self.scheme.root_mover, None)
        root = self.scheme.move_decision_tree()
        assert_not_equal(self.scheme.root_mover, None)

        assert_equal(len(root.movers), 2)
        names = ['ShootingChooser', 'RepexChooser']
        name_dict = {root.movers[i].name : i for i in range(len(root.movers))}
        for name in names:
            assert_in(name, list(name_dict.keys()))

        assert_equal(len(root.movers[name_dict['ShootingChooser']].movers), 6)
        assert_equal(len(root.movers[name_dict['RepexChooser']].movers), 4)

        new_root = self.scheme.move_decision_tree()
        assert_is(new_root, root)

        new_root = self.scheme.move_decision_tree(rebuild=True)
        assert_is_not(new_root, root)

    def test_repex_style_switching(self):
        nn_repex = NearestNeighborRepExStrategy()
        all_repex = AllSetRepExStrategy()
        default = OrganizeByMoveGroupStrategy()

        self.scheme.append([default, nn_repex])
        root = self.scheme.move_decision_tree(rebuild=True)
        assert_equal(len(self.scheme.movers['repex']), 4)

        self.scheme.append(all_repex, force=True)
        root = self.scheme.move_decision_tree(rebuild=True)
        assert_equal(len(self.scheme.movers['repex']), 6)

        self.scheme.append(nn_repex, force=True)
        root = self.scheme.move_decision_tree(rebuild=True)
        assert_equal(len(self.scheme.movers['repex']), 4)

    def test_build_balance_partners(self):
        ensA = self.scheme.network.sampling_transitions[0].ensembles[0]
        ensB = self.scheme.network.sampling_transitions[0].ensembles[1]
        hopAB = paths.EnsembleHopMover(ensemble=ensA, target_ensemble=ensB)
        hopBA = paths.EnsembleHopMover(ensemble=ensB, target_ensemble=ensA)
        self.scheme.movers['hop'] = [hopAB, hopBA]
        self.scheme.append(strategies.OrganizeByMoveGroupStrategy())
        root = self.scheme.move_decision_tree()
        self.scheme.build_balance_partners()
        assert_equal(self.scheme.balance_partners[hopAB], [hopBA])
        assert_equal(self.scheme.balance_partners[hopBA], [hopAB])

    @raises(RuntimeWarning)
    def test_build_balance_partners_premature(self):
        self.scheme.movers = {}
        self.scheme.build_balance_partners()

    @raises(RuntimeWarning)
    def test_build_balance_partners_no_partner(self):
        ensA = self.scheme.network.sampling_transitions[0].ensembles[0]
        ensB = self.scheme.network.sampling_transitions[0].ensembles[1]
        hopAB = paths.EnsembleHopMover(ensemble=ensA, target_ensemble=ensB)
        hopBA = paths.EnsembleHopMover(ensemble=ensB, target_ensemble=ensA)
        self.scheme.movers['hop'] = [hopAB]
        self.scheme.append(strategies.OrganizeByMoveGroupStrategy())
        root = self.scheme.move_decision_tree()
        self.scheme.build_balance_partners()

    @raises(RuntimeWarning)
    def test_build_balance_partners_two_partners(self):
        ensA = self.scheme.network.sampling_transitions[0].ensembles[0]
        ensB = self.scheme.network.sampling_transitions[0].ensembles[1]
        hopAB = paths.EnsembleHopMover(ensemble=ensA, target_ensemble=ensB)
        hopAB2 = paths.EnsembleHopMover(ensemble=ensA, target_ensemble=ensB)
        hopBA = paths.EnsembleHopMover(ensemble=ensB, target_ensemble=ensA)
        self.scheme.movers['hop'] = [hopAB, hopBA, hopAB2]
        self.scheme.append(strategies.OrganizeByMoveGroupStrategy())
        root = self.scheme.move_decision_tree()
        self.scheme.build_balance_partners()

    def test_sanity_check_sane(self):
        self.scheme.append([NearestNeighborRepExStrategy(),
                            OneWayShootingStrategy(),
                            OrganizeByMoveGroupStrategy()])
        root = self.scheme.move_decision_tree()
        self.scheme.sanity_check()

    @raises(AssertionError)
    def test_sanity_check_unused_sampling(self):
        ensemble_subset = self.scheme.network.sampling_transitions[0].ensembles
        self.scheme.append([
            OneWayShootingStrategy(ensembles=ensemble_subset),
            OrganizeByMoveGroupStrategy()
        ])
        root = self.scheme.move_decision_tree()
        self.scheme.sanity_check()

    @raises(AssertionError)
    def test_sanity_check_choice_prob_fails(self):
        self.scheme.append([NearestNeighborRepExStrategy(),
                            OneWayShootingStrategy(),
                            OrganizeByMoveGroupStrategy()])
        root = self.scheme.move_decision_tree()
        key0 = list(self.scheme.choice_probability.keys())[0]
        self.scheme.choice_probability[key0] = 0.0
        self.scheme.sanity_check()

    @raises(AssertionError)
    def test_sanity_check_duplicated_movers(self):
        ensemble_subset = self.scheme.network.sampling_transitions[0].ensembles
        self.scheme.append([
            OneWayShootingStrategy(),
            OrganizeByMoveGroupStrategy()
        ])
        root = self.scheme.move_decision_tree()
        self.scheme.movers['foo'] = [self.scheme.movers['shooting'][0]]
        self.scheme.sanity_check()

    @raises(TypeError)
    def test_select_movers_no_choice_probability(self):
        self.scheme.append([OneWayShootingStrategy(),
                            OrganizeByMoveGroupStrategy()])
        movers = self.scheme._select_movers('shooting')

    def test_select_movers(self):
        self.scheme.append([
            OneWayShootingStrategy(), 
            NearestNeighborRepExStrategy(),
            OrganizeByMoveGroupStrategy()
        ])
        root = self.scheme.move_decision_tree()
        some_shooters = self.scheme.movers['shooting'][0:2]

        movers = self.scheme._select_movers('shooting')
        assert_equal(movers, self.scheme.movers['shooting'])

        movers = self.scheme._select_movers(some_shooters)
        assert_equal(movers, some_shooters)

        movers = self.scheme._select_movers(some_shooters[0])
        assert_equal(movers, [some_shooters[0]])

    def test_n_trials_for_steps(self):
        self.scheme.append([
            OneWayShootingStrategy(), 
            NearestNeighborRepExStrategy(),
            OrganizeByMoveGroupStrategy()
        ])
        # we should have 6 shooters and 4 repex movers, but default strategy
        # means that have the probability of selecting repex; so we get
        # a shooting move 75% of the time
        root = self.scheme.move_decision_tree()
        root = self.scheme.move_decision_tree()
        some_shooters = self.scheme.movers['shooting'][0:2]

        assert_almost_equal(
            self.scheme.n_trials_for_steps('shooting', 100), 75.0
        )

        assert_almost_equal(
            self.scheme.n_trials_for_steps(some_shooters, 100), 25.0
        )

        assert_almost_equal(
            self.scheme.n_trials_for_steps(some_shooters[0], 100), 12.5
        )


    def test_n_steps_for_trials(self):
        self.scheme.append([
            OneWayShootingStrategy(),
            NearestNeighborRepExStrategy(),
            OrganizeByMoveGroupStrategy()
        ])
        # we should have 6 shooters and 4 repex movers, but default strategy
        # means that have the probability of selecting repex; so we get
        # a shooting move 75% of the time
        root = self.scheme.move_decision_tree()
        some_shooters = self.scheme.movers['shooting'][0:2]

        assert_almost_equal(
            self.scheme.n_steps_for_trials('shooting', 100), old_div(400.0,3.0)
        )

        assert_almost_equal(
            self.scheme.n_steps_for_trials(some_shooters, 100), 400.0
        )

        assert_almost_equal(
            self.scheme.n_steps_for_trials(some_shooters[0], 100), 800.0
        )



class TestDefaultScheme(object):
    def setup(self):
        paths.InterfaceSet._reset()
        cvA = paths.FunctionCV(name="xA", f=lambda s : s.xyz[0][0])
        cvB = paths.FunctionCV(name="xB", f=lambda s : -s.xyz[0][0])
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
        self.no_ms_outer = paths.MSTISNetwork(
            [(self.stateA, interfacesA), (self.stateB, interfacesB)]
        )

    def test_default_scheme(self):
        scheme = DefaultScheme(self.network)
        root = scheme.move_decision_tree()
        chooser_type_dict = {
            'ShootingChooser' : paths.OneWayShootingMover,
            'PathreversalChooser' : paths.PathReversalMover,
            'RepexChooser' : paths.ReplicaExchangeMover,
            'MinusChooser' : paths.MinusMover,
            'Ms_outer_shootingChooser' : paths.OneWayShootingMover
        }
        names = list(chooser_type_dict.keys())
        assert_equal(len(root.movers), len(names))

        name_dict = {root.movers[i].name : i for i in range(len(root.movers))}
        for name in names:
            assert_in(name, list(name_dict.keys()))

        n_normal_repex = 4
        n_msouter_repex = 2
        n_repex = n_normal_repex + n_msouter_repex

        assert_equal(
            len(root.movers[name_dict['ShootingChooser']].movers), 6
        )
        assert_equal(
            len(root.movers[name_dict['PathreversalChooser']].movers), 7
        )
        assert_equal(
            len(root.movers[name_dict['RepexChooser']].movers), n_repex
        )
        assert_equal(
            len(root.movers[name_dict['MinusChooser']].movers), 2
        )
        assert_equal(
            len(root.movers[name_dict['Ms_outer_shootingChooser']].movers), 1
        )

        for choosername in names:
            for mover in root.movers[name_dict[choosername]].movers:
                assert_equal(type(mover), chooser_type_dict[choosername])

    def test_default_scheme_no_ms_outer(self):
        scheme = DefaultScheme(self.no_ms_outer)
        root = scheme.move_decision_tree()
        chooser_type_dict = {
            'ShootingChooser' : paths.OneWayShootingMover,
            'PathreversalChooser' : paths.PathReversalMover,
            'RepexChooser' : paths.ReplicaExchangeMover,
            'MinusChooser' : paths.MinusMover
        }
        names = list(chooser_type_dict.keys())
        assert_equal(len(root.movers), len(names))

        name_dict = {root.movers[i].name : i for i in range(len(root.movers))}
        for name in names:
            assert_in(name, list(name_dict.keys()))

        n_normal_repex = 4
        n_msouter_repex = 0
        n_repex = n_normal_repex + n_msouter_repex

        assert_equal(
            len(root.movers[name_dict['ShootingChooser']].movers), 6
        )
        assert_equal(
            len(root.movers[name_dict['PathreversalChooser']].movers), 6
        )
        assert_equal(
            len(root.movers[name_dict['RepexChooser']].movers), n_repex
        )
        assert_equal(
            len(root.movers[name_dict['MinusChooser']].movers), 2
        )

    def test_default_sanity(self):
        scheme = DefaultScheme(self.network)
        root = scheme.move_decision_tree()
        scheme.sanity_check()

    def test_default_hidden_ensembles(self):
        scheme = DefaultScheme(self.network)
        root = scheme.move_decision_tree()
        hidden = scheme.find_hidden_ensembles()
        assert_equal(len(hidden), 2)

    def test_default_unused_ensembles(self):
        scheme = DefaultScheme(self.network)
        root = scheme.move_decision_tree()
        unused = scheme.find_unused_ensembles()
        assert_equal(len(unused), 0) # will change when minus/msouter

    def test_default_balance_partners(self):
        scheme = DefaultScheme(self.network)
        root = scheme.move_decision_tree()
        scheme.build_balance_partners()
        # by default, every mover is its own balance partner
        for group in list(scheme.movers.values()):
            for mover in group:
                assert_equal(scheme.balance_partners[mover], [mover])

    def test_default_choice_probability(self):
        scheme = DefaultScheme(self.network)
        root = scheme.move_decision_tree()
        default_group_weights = {
            'shooting' : 1.0,
            'repex' : 0.5,
            'pathreversal' : 0.5,
            'minus' : 0.2,
            'ms_outer_shooting' : 1.0
        }

        assert_almost_equal(sum(scheme.choice_probability.values()), 1.0)

        tot_norm = sum([default_group_weights[group] 
                        for group in scheme.movers])

        prob_shoot0 = scheme.choice_probability[scheme.movers['shooting'][0]]
        n_shooting = len(scheme.movers['shooting'])
        for group in default_group_weights:
            scale = default_group_weights[group]
            n_group = len(scheme.movers[group])
            expected_prob = default_group_weights[group]*prob_shoot0
            for mover in scheme.movers[group]:
                test_prob = scheme.choice_probability[mover]
                assert_almost_equal(expected_prob, test_prob)

    def test_initial_conditions_from_trajectory(self):
        scheme = DefaultScheme(self.network)
        # root = scheme.move_decision_tree()
        assert_equal(len(scheme.list_initial_ensembles()), 9)

        traj1 = make_1d_traj([-0.6, -0.2, -0.6])
        traj2 = make_1d_traj([-0.6, -0.2, -0.05, -0.4, -0.6])
        traj3 = make_1d_traj([-0.6, -0.2, 0.2, 0.6])

        all_trajs = [traj1, traj2, traj3]

        traj1r = traj1.reversed
        traj2r = traj2.reversed
        traj3r = traj3.reversed

        def assert_init_cond(sample_set, ensembles, expected):
            # helper to check the results. Expected is in the form
            # of a list of resulting trajectories
            # ens is the list of ensembles to be tested in order

            sample_set.sanity_check()

            assert_equal(len(sample_set), len(expected))

            for ensemble, traj in zip(ensembles, expected):
                # print ensemble.name, sample_set[ensemble].trajectory.xyz[:,0,0], traj.xyz[:, 0,0],
                # print hex(id(traj)), hex(id(sample_set[ensemble].trajectory.xyz[:,0,0]))
                assert_equal(sample_set[ensemble].trajectory, traj)

        transAB = transBA = None
        for trans in self.network.sampling_transitions:
            if trans.stateA == self.stateA and trans.stateB == self.stateB:
                transAB = trans
            elif trans.stateA == self.stateB and trans.stateB == self.stateA:
                transBA = trans
            else:
                raise RuntimeWarning("That's a weird transition!")

        ms_outer_ens = list(self.network.special_ensembles['ms_outer'].keys())[0]

        ensembles = transAB.ensembles + [ms_outer_ens] + transBA.ensembles

        init_cond = scheme.initial_conditions_from_trajectories(
            trajectories=all_trajs,
            preconditions=[],
            reuse_strategy='all',
            strategies=['get']
        )

        assert_init_cond(
            init_cond,
            ensembles[:4],
            [traj1, traj1, traj2, traj3]
        )

        init_cond = scheme.initial_conditions_from_trajectories(
            trajectories=all_trajs,
            preconditions=['mirror'],
            reuse_strategy='all',
            strategies=['get']
        )

        assert_init_cond(
            init_cond, ensembles,
            [traj1, traj1, traj2] + [traj3] + [traj3r] * 3
        )

        init_cond = scheme.initial_conditions_from_trajectories(
            trajectories=all_trajs,
            preconditions=['mirror', 'sort-shortest'],
            strategies=['get'],
            reuse_strategy='all'
        )
        assert_init_cond(
            init_cond, ensembles,
            [traj1, traj1, traj3] + [traj3] + [traj3r] * 3
        )

        init_cond = scheme.initial_conditions_from_trajectories(
            trajectories=all_trajs,
            preconditions=[],
            reuse_strategy='avoid',
            strategies=['get']
        )

        assert_init_cond(
            init_cond,
            ensembles[:4],
            [traj1, traj2, traj3, traj3]
        )

        init_cond = scheme.initial_conditions_from_trajectories(
            trajectories=all_trajs,
            preconditions=['mirror'],
            reuse_strategy='avoid',
            strategies=['get']
        )

        init_cond.sanity_check()
        assert_equal(len(init_cond), 7)

        for ensemble, traj in zip(ensembles[:3], [traj1, traj2, traj3]):
            assert_equal(init_cond[ensemble].trajectory, traj)
        for ensemble, traj in zip(ensembles[4:], [traj3r] * 3):
            assert_equal(init_cond[ensemble].trajectory, traj)

        # because of the way the scheme ensembles are creating involving a
        # set, the order in which the ensemble are created changes.
        # in some cases traj3 is used and hence avoided in the outer
        # in some cases traj3r, but both are fine.
        try:
            assert_equal(init_cond[ensembles[3]].trajectory, traj3)
        except AssertionError:
            assert_equal(init_cond[ensembles[3]].trajectory, traj3r)

        init_cond = scheme.initial_conditions_from_trajectories(
            trajectories=all_trajs,
            preconditions=['mirror', 'sort-shortest'],
            reuse_strategy='avoid',
            strategies=['get']
        )
        init_cond.sanity_check()
        assert_equal(len(init_cond), 7)

        for ensemble, traj in zip(ensembles[:3], [traj1, traj1r, traj3]):
            assert_equal(init_cond[ensemble].trajectory, traj)
        for ensemble, traj in zip(ensembles[4:], [traj3r] * 3):
            assert_equal(init_cond[ensemble].trajectory, traj)

        try:
            assert_equal(init_cond[ensembles[3]].trajectory, traj3)
        except AssertionError:
            assert_equal(init_cond[ensembles[3]].trajectory, traj3r)

        # this one avoids reversed copies
        init_cond = scheme.initial_conditions_from_trajectories(
            trajectories=[traj1],
            preconditions=[],
            strategies=['get']
        )
        assert_init_cond(
            init_cond,
            ensembles[:2],
            [traj1, traj1]
        )
        init_cond = scheme.initial_conditions_from_trajectories(
            trajectories=traj2,
            sample_set=init_cond,
            preconditions=[],
            strategies=['get']
        )

        assert_init_cond(
            init_cond,
            ensembles[:3],
            [traj1, traj1, traj2]
        )

        init_cond = scheme.initial_conditions_from_trajectories(
            trajectories=[traj3],
            preconditions=[],
            strategies=['get']
        )

        assert_init_cond(
            init_cond,
            ensembles[:4],
            [traj3] * 4
        )

        init_cond = scheme.initial_conditions_from_trajectories(
            trajectories=[traj3],
            preconditions=['mirror'],
            strategies=['get']
        )
        assert_init_cond(
            init_cond,
            ensembles,
            [traj3] * 4 + [traj3r] * 3
        )

    def test_check_initial_conditions(self):
        scheme = DefaultScheme(self.network)
        traj3 = make_1d_traj([-0.6, -0.2, 0.2, 0.6])
        # cheating a bit, since we know what this gives
        init_cond = scheme.initial_conditions_from_trajectories(traj3)
        assert_equal(len(init_cond), 7)
        assert_equal(len(scheme.list_initial_ensembles()), 9)
        (missing, extra) = scheme.check_initial_conditions(init_cond)
        assert_equal(len(missing), 2)
        assert_equal(len(extra), 0)
        for ens in list(self.network.special_ensembles['minus'].keys()):
            assert_in([ens], missing)
        init_cond.append_as_new_replica(
            paths.Sample(trajectory=traj3,
                         ensemble=paths.LengthEnsemble(4),
                         replica=None)
        )
        (missing, extra) = scheme.check_initial_conditions(init_cond)
        assert_equal(len(missing), 2)
        assert_equal(len(extra), 1)

    @raises(AssertionError)
    def test_assert_initial_conditions(self):
        scheme = DefaultScheme(self.network)
        traj3 = make_1d_traj([-0.6, -0.2, 0.2, 0.6])
        init_cond = scheme.initial_conditions_from_trajectories(traj3)
        init_cond.append_as_new_replica(
            paths.Sample(trajectory=traj3,
                         ensemble=paths.LengthEnsemble(4),
                         replica=None)
        )
        scheme.assert_initial_conditions(init_cond)

    def test_initial_conditions_report(self):
        scheme = DefaultScheme(self.network)
        traj3 = make_1d_traj([-0.6, -0.2, 0.2, 0.6])
        init_cond = scheme.initial_conditions_from_trajectories(traj3)
        init_cond.append_as_new_replica(
            paths.Sample(trajectory=traj3,
                         ensemble=paths.LengthEnsemble(4),
                         replica=None)
        )
        start = "Missing ensembles:\n"
        missing_A = "*  [Out A minus]\n"
        missing_B = "*  [Out B minus]\n"
        finish = "Extra ensembles:\n*  [LengthEnsemble]\n"
        expected_AB = start + missing_A + missing_B + finish
        expected_BA = start + missing_B + missing_A + finish
        result = scheme.initial_conditions_report(init_cond)
        try:
            assert_equal(result, expected_AB)
        except AssertionError:
            assert_equal(result, expected_BA)


class TestLockedMoveScheme(object):
    def setup(self):
        paths.InterfaceSet._reset()
        cvA = paths.FunctionCV(name="xA", f=lambda s : s.xyz[0][0])
        cvB = paths.FunctionCV(name="xB", f=lambda s : -s.xyz[0][0])
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
        self.basic_scheme = DefaultScheme(self.network)
        self.root_mover = self.basic_scheme.move_decision_tree()

    def test_initialization(self):
        scheme = LockedMoveScheme(self.root_mover, self.network)
        assert_equal(scheme.network, self.network)
        assert_equal(scheme.move_decision_tree(), self.root_mover)

    def test_build_move_decision_tree(self):
        scheme = LockedMoveScheme(self.root_mover, self.network)
        scheme.move_decision_tree(rebuild=True)
        assert_equal(scheme.move_decision_tree(), self.root_mover)

    @raises(TypeError)
    def test_append(self):
        scheme = LockedMoveScheme(self.root_mover, self.network)
        scheme.append(AllSetRepExStrategy())

    @raises(TypeError)
    def test_apply_strategy(self):
        scheme = LockedMoveScheme(self.root_mover, self.network)
        strategy = AllSetRepExStrategy()
        scheme.apply_strategy(strategy)

    @raises(AttributeError)
    def test_choice_probability_fail(self):
        scheme = LockedMoveScheme(self.root_mover, self.network)
        vals = scheme.choice_probability

    def test_choice_probability_works(self):
        scheme = LockedMoveScheme(self.root_mover, self.network)
        scheme.choice_probability = self.basic_scheme.choice_probability
        vals = scheme.choice_probability

    @raises(AttributeError)
    def test_movers_fail(self):
        scheme = LockedMoveScheme(self.root_mover, self.network)
        vals = scheme.movers

    def test_movers_works(self):
        scheme = LockedMoveScheme(self.root_mover, self.network)
        scheme.movers = self.basic_scheme.movers
        vals = scheme.movers


class TestOneWayShootingMoveScheme(object):
    def setup(self):
        paths.InterfaceSet._reset()
        cvA = paths.FunctionCV(name="xA", f=lambda s : s.xyz[0][0])
        cvB = paths.FunctionCV(name="xB", f=lambda s : -s.xyz[0][0])
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

    def test_scheme(self):
        scheme = OneWayShootingMoveScheme(self.network)
        root = scheme.move_decision_tree()
        assert_equal(len(scheme.movers), 1)
        assert_equal(len(root.movers), 1)

    def test_sanity(self):
        scheme = OneWayShootingMoveScheme(self.network)
        root = scheme.move_decision_tree()
        scheme.sanity_check()

    def test_unused_ensembles(self):
        scheme = OneWayShootingMoveScheme(self.network)
        root = scheme.move_decision_tree()
        unused = scheme.find_unused_ensembles()
        specials = self.network.special_ensembles
        expected_unused = sum([list(specials[special_type].keys()) 
                               for special_type in specials], [])
        assert_equal(set(expected_unused), set(unused))

    def test_check_initial_conditions(self):
        scheme = OneWayShootingMoveScheme(self.network)
        traj3 = make_1d_traj([-0.6, -0.2, 0.2, 0.6])
        init_cond = scheme.initial_conditions_from_trajectories(traj3)
        assert_equal(len(scheme.list_initial_ensembles()), 6)
        assert_equal(len(init_cond), 6)
        scheme.assert_initial_conditions(init_cond)
        assert_equal(scheme.initial_conditions_report(init_cond),
                     "No missing ensembles.\nNo extra ensembles.\n")
