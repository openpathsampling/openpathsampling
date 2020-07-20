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

import math

import openpathsampling as paths
from openpathsampling.high_level.move_scheme import *
from openpathsampling.high_level.move_strategy import (
    levels,
    MoveStrategy, OneWayShootingStrategy, NearestNeighborRepExStrategy,
    OrganizeByMoveGroupStrategy, AllSetRepExStrategy
)

import openpathsampling.high_level.move_strategy as strategies

import pytest

import logging
logging.getLogger('openpathsampling.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.ensemble').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.netcdfplus').setLevel(logging.CRITICAL)

def _make_acceptance_mock_step(mccycle, accepted, path_sim_mover, move_type,
                               mover_sig, submover_num=None):
    root_mover = path_sim_mover.mover
    chooser_names = {m.name[:-7].lower(): m for m in root_mover.movers}
    chooser = chooser_names[move_type]
    sig_to_mover = {frozenset(m.ensemble_signature[0]): m
                    for m in chooser.movers}
    # group_mover = chooser.movers[mover_num]
    group_mover = sig_to_mover[frozenset(mover_sig)]
    Change = {True: paths.AcceptedSampleMoveChange,
              False: paths.RejectedSampleMoveChange}[accepted]
    # foo here is because we need non-empty samples to show that we're
    # actually accepted or not (WHY?!?!?)
    if submover_num is not None:
        submover = group_mover.movers[submover_num]
        submover_change = Change(samples=['foo'], mover=submover)
        group_mover_change = paths.RandomChoiceMoveChange(
            subchange=submover_change,
            mover=group_mover
        )
    else:
        submover_change = None
        group_mover_change = Change(samples=['foo'], mover=group_mover)

    chooser_change = paths.RandomChoiceMoveChange(
        subchange=group_mover_change,
        mover=chooser
    )
    root_mover_change = paths.RandomChoiceMoveChange(
        subchange=chooser_change,
        mover=root_mover
    )
    path_sim_change = paths.PathSimulatorMoveChange(
        subchange=root_mover_change,
        mover=path_sim_mover
    )
    step = paths.MCStep(
        mccycle=mccycle,
        active=paths.SampleSet([]),
        change=path_sim_change
    )
    return step

def _make_null_mover_step(mccycle, path_sim_mover, null_mover):
    empty_sample_set = paths.SampleSet([])
    change = paths.PathSimulatorMoveChange(
        mover=path_sim_mover,
        subchange=null_mover.move(empty_sample_set)
    )
    step = paths.MCStep(
        mccycle=mccycle,
        active=empty_sample_set,
        change=change
    )
    return step


class TestMoveAcceptanceAnalysis(object):
    def setup(self):
        self.HAS_TQDM = paths.progress.HAS_TQDM
        paths.progress.HAS_TQDM = False
        paths.InterfaceSet._reset()
        cvA = paths.FunctionCV(name="xA", f=lambda s : s.xyz[0][0])
        cvB = paths.FunctionCV(name="xB", f=lambda s : -s.xyz[0][0])
        state_A = paths.CVDefinedVolume(cvA, float("-inf"), -0.5).named("A")
        state_B = paths.CVDefinedVolume(cvB, float("-inf"), -0.5).named("B")
        interfaces_A = paths.VolumeInterfaceSet(cvA, float("-inf"),
                                               [-0.5, -0.3])
        network = paths.MISTISNetwork([(state_A, interfaces_A, state_B)])
        self.scheme = MoveScheme(network)
        self.scheme.append(OneWayShootingStrategy())
        self.scheme.append(NearestNeighborRepExStrategy())
        self.scheme.append(OrganizeByMoveGroupStrategy())

        root_mover = self.scheme.move_decision_tree()
        path_sim_mover = paths.PathSimulatorMover(root_mover, None)
        null_mover = paths.IdentityPathMover(counts_as_trial=False)

        ens_0 = network.sampling_ensembles[0]
        ens_1 = network.sampling_ensembles[1]

        # acc repex ens1-2
        # acc   fwd ens1
        # acc  bkwd ens2
        # rej  bkwd ens1
        # rej repex ens1-2
        step_info = [
            (1, True, path_sim_mover, 'repex', [ens_0, ens_1], None),
            (2, True, path_sim_mover, 'shooting', [ens_0], 0),
            (3, True, path_sim_mover, 'shooting', [ens_1], 1),
            (4, False, path_sim_mover, 'shooting', [ens_0], 1),
            (5, False, path_sim_mover, 'repex', [ens_0, ens_1], None)
        ]
        self.steps = [_make_acceptance_mock_step(*info)
                      for info in step_info]

        self.null_mover_6 = _make_null_mover_step(6, path_sim_mover,
                                                  null_mover)
        self.null_mover_change_key = [(None, str([path_sim_mover, [None]]))]

        acceptance_empty = MoveAcceptanceAnalysis(self.scheme)

        acceptance = MoveAcceptanceAnalysis(self.scheme)
        acceptance.add_steps(self.steps)

        acceptance_null = MoveAcceptanceAnalysis(self.scheme)
        acceptance_null.add_steps(self.steps + [self.null_mover_6])

        self.analysis = {'empty': acceptance_empty,
                         'normal': acceptance,
                         'with_null': acceptance_null}

    def teardown(self):
        paths.progress.HAS_TQDM = self.HAS_TQDM

    @pytest.mark.parametrize('step_num', [0, 1, 2, 3, 4])
    def test_calculate_step_acceptance(self, step_num):
        accepted = [1] if step_num in [0, 1, 2] else [0]
        analysis = MoveAcceptanceAnalysis(self.scheme)
        analysis._calculate_step_acceptance(self.steps[step_num])
        assert len(analysis._trials) == len(analysis._accepted)
        len_trials = len(analysis._trials)
        assert list(analysis._trials.values()) == [1] * len_trials
        assert list(analysis._accepted.values()) == accepted * len_trials

    def test_add_steps(self):
        # also tests n_total_trials
        acceptance = MoveAcceptanceAnalysis(self.scheme)
        assert acceptance._n_steps == 0
        assert acceptance.n_total_trials == 0
        acceptance.add_steps(self.steps)
        assert acceptance._n_steps == 5
        assert acceptance.n_total_trials == 5
        acceptance.add_steps([self.null_mover_6])
        assert acceptance._n_steps == 6
        assert acceptance.n_total_trials == 5

    @pytest.mark.parametrize('simulation', ['empty', 'normal', 'with_null'])
    def test_no_move_keys(self, simulation):
        analysis = self.analysis[simulation]
        expected = {'empty': [],
                    'normal': [],
                    'with_null': self.null_mover_change_key}[simulation]
        assert analysis.no_move_keys == expected

    def test_select_movers_none(self):
        analysis = self.analysis['normal']  # doesn't matter which
        select_movers = analysis._select_movers
        scheme_movers = {k: self.scheme.movers[k]
                         for k in ['shooting', 'repex']}
        assert select_movers(None) == scheme_movers

    @pytest.mark.parametrize('group_name', ['shooting', 'repex'])
    def test_select_movers_groupname(self, group_name):
        analysis = self.analysis['normal']  # doesn't matter which
        select_movers = analysis._select_movers
        expected = {mover: [mover]
                    for mover in self.scheme.movers[group_name]}

        assert select_movers(group_name) == expected

    @pytest.mark.parametrize('group_name', ['shooting', 'repex'])
    def test_select_movers_mover(self, group_name):
        analysis = self.analysis['normal']  # doesn't matter which
        select_movers = analysis._select_movers
        input_movers = self.scheme.movers[group_name]
        for mover in input_movers:
            try:
                extra_movers = mover.movers
            except AttributeError:
                extra_movers = []

            expected = {m: [m] for m in [mover] + extra_movers}
            assert select_movers(mover) == expected

    @pytest.mark.parametrize('simulation', ['empty', 'normal', 'with_null'])
    def test_summary_data_none(self, simulation):
        results = self.analysis[simulation].summary_data(None)
        expected_results_empty = {
            'shooting': {'move_name': 'shooting',
                         'n_accepted': 0,
                         'n_trials': 0,
                         'expected_frequency': 0.8},
            'repex': {'move_name': 'repex',
                      'n_accepted': 0,
                      'n_trials': 0,
                      'expected_frequency': 0.2}
        }
        expected_results_non_empty = {
            'shooting': {'move_name': 'shooting',
                         'n_accepted': 2,
                         'n_trials': 3,
                         'expected_frequency': 0.8},
            'repex': {'move_name': 'repex',
                      'n_accepted': 1,
                      'n_trials': 2,
                      'expected_frequency': 0.2}
        }
        expected = {'empty': expected_results_empty,
                    'normal': expected_results_non_empty,
                    'with_null': expected_results_non_empty}[simulation]

        for result in results:
            assert result._asdict() == expected[result.move_name]

    @pytest.mark.parametrize('group_name', ['shooting', 'repex'])
    @pytest.mark.parametrize('simulation', ['empty', 'normal', 'with_null'])
    def test_summary_data_groupname(self, group_name, simulation):
        results = self.analysis[simulation].summary_data(group_name)
        scheme = self.scheme

        analysis = self.analysis[simulation]
        for i, mover in enumerate(scheme.movers['shooting']):
            print(i, mover, [v for k, v in  analysis._trials.items()
                             if k[0] == mover])
        expected_results_empty = {
            'shooting': [{'move_name': scheme.movers['shooting'][0],
                          'expected_frequency': 0.4,
                          'n_accepted': 0,
                          'n_trials': 0},
                         {'move_name': scheme.movers['shooting'][1],
                          'expected_frequency': 0.4,
                          'n_accepted': 0,
                          'n_trials': 0}],
            'repex': [{'move_name': scheme.movers['repex'][0],
                       'expected_frequency': 0.2,
                       'n_accepted': 0,
                       'n_trials': 0}]
        }
        expected_results_non_empty = {
            'shooting': [{'move_name': scheme.movers['shooting'][0],
                          'expected_frequency': 0.4,
                          'n_accepted': 1,
                          'n_trials': 2},
                         {'move_name': scheme.movers['shooting'][1],
                          'expected_frequency': 0.4,
                          'n_accepted': 1,
                          'n_trials': 1}],
            'repex': [{'move_name': scheme.movers['repex'][0],
                       'expected_frequency': 0.2,
                       'n_accepted': 1,
                       'n_trials': 2}]
        }
        expected_list = {'empty': expected_results_empty,
                         'normal': expected_results_non_empty,
                         'with_null': expected_results_non_empty
                        }[simulation][group_name]
        expected = {res['move_name']: res for res in expected_list}

        for result in results:
            assert result._asdict() == expected[result.move_name]

    @pytest.mark.parametrize('simulation', ['empty', 'normal', 'with_null'])
    @pytest.mark.parametrize('mover_ensemble', [0, 1])
    def test_summary_data_mover(self, simulation, mover_ensemble):
        mover = self.scheme.movers['shooting'][mover_ensemble]
        results = self.analysis[simulation].summary_data(mover)

        expected_list = [
            {'move_name': mover,
             # 'expected_frequency': 0.4,
             'n_accepted': 0,
             'n_trials': 0},
            {'move_name': mover.movers[0],
             # 'expected_frequency': float('nan'),
             'n_accepted': 0,
             'n_trials': 0},
            {'move_name': mover.movers[1],
             # 'expected_frequency': float('nan'),
             'n_accepted': 0,
             'n_trials': 0}
        ]
        updates = {
            0: {mover: {'n_accepted': 1, 'n_trials': 2},
                mover.movers[0]: {'n_accepted': 1, 'n_trials': 1},
                mover.movers[1]: {'n_accepted': 0, 'n_trials': 1}},
            1: {mover: {'n_accepted': 1, 'n_trials': 1},
                mover.movers[0]: {'n_accepted': 0, 'n_trials': 0},
                mover.movers[1]: {'n_accepted': 1, 'n_trials': 1}}
        }

        expected = {elem['move_name']: elem for elem in expected_list}
        if simulation in ['normal', 'with_null']:
            update = updates[mover_ensemble]
            for m, dct in expected.items():
                dct.update(update[m])

        for result in results:
            # trickiness here because 'nan' != 'nan'
            result_dict = result._asdict()
            result_freq = result_dict.pop('expected_frequency')
            if result.move_name == mover:
                assert result_freq == 0.4
            else:
                assert math.isnan(result_freq)

            assert result_dict == expected[result.move_name]

    @pytest.mark.parametrize('simulation', ['normal', 'with_null'])
    def test_line_as_text(self, simulation):
        line = MoveAcceptanceAnalysisLine(
            move_name='shooting',
            n_accepted=2,
            n_trials=3,
            expected_frequency=0.8
        )
        expected = ("shooting ran 60.000% (expected 80.00%) of the "
                    + "cycles with acceptance 2/3 (66.67%)\n")
        result = self.analysis[simulation]._line_as_text(line)
        assert result == expected

    def test_line_as_text_nan_acceptance(self):
        line = MoveAcceptanceAnalysisLine(
            move_name='path_reversal',
            n_accepted=0,
            n_trials=0,
            expected_frequency=float('nan')
        )
        expected = ("path_reversal ran 0.000% (expected nan%) of the "
                    + "cycles with acceptance 0/0 (nan%)\n")
        result = self.analysis['normal']._line_as_text(line)
        assert result == expected

    def test_line_as_text_mover_as_name(self):
        mover = self.scheme.movers['shooting'][0]
        line = MoveAcceptanceAnalysisLine(
            move_name=mover,
            n_accepted=1,
            n_trials=2,
            expected_frequency=0.4
        )
        expected = (str(mover) + " ran 40.000% (expected 40.00%) "
                    + "of the cycles with acceptance 1/2 (50.00%)\n")
        result = self.analysis['normal']._line_as_text(line)
        assert result == expected

    @pytest.mark.parametrize('simulation', ['normal', 'with_null'])
    def test_format_as_text(self, simulation):
        analysis = self.analysis[simulation]
        summary_data = analysis.summary_data(None)
        text_lines = {
            'shooting': ("shooting ran 60.000% (expected 80.00%) of the "
                       + "cycles with acceptance 2/3 (66.67%)\n"),
            'repex': ("repex ran 40.000% (expected 20.00%) of the "
                      + "cycles with acceptance 1/2 (50.00%)\n")
        }
        expected = "".join([text_lines[line.move_name]
                            for line in summary_data])

        if simulation == 'with_null':
            expected = ("Null moves for 1 cycles. Excluding null moves:\n"
                        + expected)

        assert analysis.format_as_text(summary_data) == expected


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
