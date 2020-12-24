import itertools
import random
import pytest
from nose.tools import assert_equal, assert_almost_equal, raises
from .test_helpers import (make_1d_traj, MoverWithSignature, RandomMDEngine,
                           assert_frame_equal, assert_items_equal)

from openpathsampling.analysis.tis import *
from openpathsampling.analysis.tis.core import steps_to_weighted_trajectories
from openpathsampling.analysis.tis.flux import default_flux_sort
import openpathsampling as paths

import pandas as pd
import pandas.testing as pdt

import logging
logging.getLogger('openpathsampling.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.ensemble').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.netcdfplus').setLevel(logging.CRITICAL)


def make_tis_traj_fixed_steps(n_steps, step_size=0.1, reverse=False):
    if reverse:
        sign = -1
        delta = 1.0
    else:
        sign = 1
        delta = 0.0
    rising = [delta + sign * (-0.5 + i) * step_size
              for i in range(n_steps + 2)]
    falling = list(reversed(rising))[1:]
    return make_1d_traj(rising + falling)


class TestMultiEnsembleSamplingAnalyzer(object):
    # this only has to check the error generation; everything else gets
    # covered by tests of subclasses
    @raises(RuntimeError)
    def test_no_ensembles(self):
        histogrammer = MultiEnsembleSamplingAnalyzer()
        histogrammer.calculate([])


class TISAnalysisTester(object):
    # abstract class to give the same setup to all the test functions

    def _make_fake_steps(self, sample_sets, mover):
        steps = []
        for (mccycle, sample_set) in enumerate(sample_sets):
            change = paths.AcceptedSampleMoveChange(
                samples=sample_set.samples,
                mover=mover,
                details=None,
                input_samples=None
            )
            step = paths.MCStep(mccycle=mccycle,
                                active=sample_set,
                                change=change)
            steps.append(step)
        return steps

    def _make_fake_sampling_sets(self, network):
        ensembles_AB = self.sampling_ensembles_for_transition(
            network, self.state_A, self.state_B
        )
        ensembles_BA = self.sampling_ensembles_for_transition(
            network, self.state_B, self.state_A
        )

        all_ensembles = ensembles_AB + ensembles_BA

        # This encodes how the SampleSets are at each time step. This is the
        # trajectory number (from trajs_AB/trajs_BA) for each ensemble
        # (index order of sampling_AB.ensembles + sampling_BA.ensembles)
        descriptions = [
            [2, 3, 3, 2, 3, 3],
            [1, 2, 3, 1, 2, 3],
            [0, 1, 2, 0, 1, 2],
            [0, 1, 2, 0, 1, 2]
        ]

        # here's the fancy fake data
        sample_sets = []
        for descr in descriptions:
            set_trajectories = ([self.trajs_AB[d] for d in descr[:3]]
                                + [self.trajs_BA[d] for d in descr[3:]])
            sample_set = paths.SampleSet([
                paths.Sample(trajectory=traj,
                             ensemble=ens,
                             replica=rep)
                for (traj, ens, rep) in zip(set_trajectories, all_ensembles,
                                            range(len(all_ensembles)))
            ])
            sample_sets.append(sample_set)
        return sample_sets

    def sampling_ensembles_for_transition(self, network, state_A, state_B):
        analysis_AB = network.transitions[(state_A, state_B)]
        sampling_AB = network.analysis_to_sampling[analysis_AB][0]
        return sampling_AB.ensembles

    def setup(self):
        # set up the trajectories, ensembles, etc. for this test
        self.HAS_TQDM = paths.progress.HAS_TQDM
        paths.progress.HAS_TQDM = False  # turn of progress bars
        paths.InterfaceSet._reset()
        cv_A = paths.FunctionCV('Id', lambda s: s.xyz[0][0])
        cv_B = paths.FunctionCV('1-Id', lambda s: 1.0-s.xyz[0][0])
        self.cv_x = cv_A
        self.state_A = paths.CVDefinedVolume(cv_A,
                                             float("-inf"), 0.0).named("A")
        self.state_B = paths.CVDefinedVolume(cv_B,
                                             float("-inf"), 0.0).named("B")
        interfaces_AB = paths.VolumeInterfaceSet(cv_A, float("-inf"),
                                                 [0.0, 0.1, 0.2])
        interfaces_BA = paths.VolumeInterfaceSet(cv_B, float("-inf"),
                                                 [0.0, 0.1, 0.2])

        # trajectory that crosses each interface, one state-to-state
        self.trajs_AB = [make_tis_traj_fixed_steps(i) for i in [0, 1, 2]]
        self.trajs_AB += [make_1d_traj([(-0.5 + i) * 0.1
                                        for i in range(12)])]

        self.trajs_BA = [make_tis_traj_fixed_steps(i, reverse=True)
                         for i in [0, 1, 2]]
        self.trajs_BA += [make_1d_traj([1.0 - (-0.5 + i) * 0.1
                                        for i in range(12)])]

        # set up mistis
        self.mistis = paths.MISTISNetwork([
            (self.state_A, interfaces_AB, self.state_B),
            (self.state_B, interfaces_BA, self.state_A)
        ])
        mover_stub_mistis = MoverWithSignature(self.mistis.all_ensembles,
                                               self.mistis.all_ensembles)

        mistis_ssets = self._make_fake_sampling_sets(self.mistis)
        self.mistis_steps = self._make_fake_steps(mistis_ssets,
                                                  mover_stub_mistis)

        self.mistis_weighted_trajectories = steps_to_weighted_trajectories(
            self.mistis_steps,
            self.mistis.sampling_ensembles
        )

        # TODO: set up mstis
        self.mstis = paths.MSTISNetwork([
            (self.state_A, interfaces_AB),
            (self.state_B, interfaces_BA)
        ])
        mover_stub_mstis = MoverWithSignature(self.mstis.all_ensembles,
                                              self.mstis.all_ensembles)
        mstis_ssets = self._make_fake_sampling_sets(self.mstis)
        self.mstis_steps = self._make_fake_steps(mstis_ssets,
                                                 mover_stub_mstis)

        self.mstis_weighted_trajectories = steps_to_weighted_trajectories(
            self.mstis_steps,
            self.mstis.sampling_ensembles
        )


    def teardown(self):
        paths.progress.HAS_TQDM = self.HAS_TQDM


class TestWeightedTrajectories(TISAnalysisTester):
    def _check_network_results(self, network, weighted_trajs):
        # works for both MISTIS and MSTIS, since they use equivalent data
        ensembles_AB = self.sampling_ensembles_for_transition(
            network, self.state_A, self.state_B
        )
        ensembles_BA = self.sampling_ensembles_for_transition(
            network, self.state_B, self.state_A
        )

        # (ensemble_number, trajectory_number): count
        results = {(0, 0): 2, (0, 1): 1, (0, 2): 1, (0, 3): 0,
                   (1, 0): 0, (1, 1): 2, (1, 2): 1, (1, 3): 1,
                   (2, 0): 0, (2, 1): 0, (2, 2): 2, (2, 3): 2}

        for ((ens, traj), result) in results.items():
            assert_equal(
                weighted_trajs[ensembles_AB[ens]][self.trajs_AB[traj]],
                result
            )
            assert_equal(
                weighted_trajs[ensembles_BA[ens]][self.trajs_BA[traj]],
                result
            )

    def test_steps_to_weighted_trajectories(self):
        assert_equal(len(self.mistis_weighted_trajectories),
                     len(self.mistis.sampling_ensembles))
        self._check_network_results(self.mistis,
                                    self.mistis_weighted_trajectories)

        assert_equal(len(self.mstis_weighted_trajectories),
                     len(self.mstis.sampling_ensembles))
        self._check_network_results(self.mstis,
                                    self.mstis_weighted_trajectories)


class TestFluxToPandas(TISAnalysisTester):
    # includes tests for default_flux_sort and flux_matrix_pd
    # as a class to simplify setup of flux objects
    def setup(self):
        super(TestFluxToPandas, self).setup()
        interfaces_A = self.mstis.from_state[self.state_A].interfaces
        interfaces_B = self.mstis.from_state[self.state_B].interfaces
        pairs_A = list(itertools.product([self.state_A], interfaces_A))
        pairs_B = list(itertools.product([self.state_B], interfaces_B))
        # note that this gives the canonical order we desire
        self.all_pairs = pairs_A + pairs_B
        self.default_ordered_results = [
            ((self.state_A, interfaces_A[0]), 10.0),
            ((self.state_A, interfaces_A[1]), 5.0),
            ((self.state_A, interfaces_A[2]), 2.0),
            ((self.state_B, interfaces_B[0]), 1.0),
            ((self.state_B, interfaces_B[1]), 0.5),
            ((self.state_B, interfaces_B[2]), 0.2)
        ]
        shuffled_fluxes = self.default_ordered_results[:]
        random.shuffle(shuffled_fluxes)
        self.fluxes = {key: value for (key, value) in shuffled_fluxes}
        self.indices = [
            ("A", "-inf<Id<0.0"), ("A", "-inf<Id<0.1"), ("A", "-inf<Id<0.2"),
            ("B", "-inf<1-Id<0.0"), ("B", "-inf<1-Id<0.1"),
            ("B", "-inf<1-Id<0.2")
        ]
        values = [value for (key, value) in self.default_ordered_results]
        index = pd.MultiIndex.from_tuples(self.indices,
                                          names=["State", "Interface"])
        self.expected_series = pd.Series(values, name="Flux", index=index)

    def test_default_flux_sort(self):
        shuffled = self.all_pairs[:]
        random.shuffle(shuffled)
        sorted_result = default_flux_sort(shuffled)
        assert_items_equal(sorted_result, self.all_pairs)

    def test_flux_matrix_pd_default(self):
        series = flux_matrix_pd(self.fluxes)
        pdt.assert_series_equal(series, self.expected_series)

    def test_flux_matrix_pd_None(self):
        series = flux_matrix_pd(self.fluxes, sort_method=None)
        for idx in self.indices:
            assert_almost_equal(series[idx], self.expected_series[idx])

    @raises(KeyError)
    def test_flux_matrix_pd_unknown_str(self):
        flux_matrix_pd(self.fluxes, sort_method="foo")


class TestDictFlux(TISAnalysisTester):
    def setup(self):
        super(TestDictFlux, self).setup()
        self.innermost_interface_A = \
            self.sampling_ensembles_for_transition(self.mistis,
                                                   self.state_A,
                                                   self.state_B)[0]
        self.innermost_interface_B = \
            self.sampling_ensembles_for_transition(self.mistis,
                                                   self.state_B,
                                                   self.state_A)[0]

        self.flux_dict = {(self.state_A, self.innermost_interface_A): 1.0,
                          (self.state_B, self.innermost_interface_B): 1.0}
        self.flux_method = DictFlux(self.flux_dict)

    def test_calculate(self):
        assert_equal(self.flux_method.calculate(self.mistis_steps),
                     self.flux_dict)

    def test_from_weighted_trajectories(self):
        assert_equal(
            self.flux_method.from_weighted_trajectories(self.mistis_steps),
            self.flux_dict
        )

    def test_intermediates(self):
        assert_equal(self.flux_method.intermediates(self.mistis_steps), [])

    def test_calculate_from_intermediates(self):
        intermediates = self.flux_method.intermediates(self.mistis_steps)
        assert_equal(
            self.flux_method.calculate_from_intermediates(*intermediates),
            self.flux_dict
        )

    def test_combine_results(self):
        my_result = self.flux_method.calculate(self.mistis_steps)
        same_result = {(self.state_A, self.innermost_interface_A): 1.0,
                       (self.state_B, self.innermost_interface_B): 1.0}
        assert_equal(
            self.flux_method.combine_results(my_result, same_result),
            my_result
        )

    @raises(RuntimeError)
    def test_bad_combine_results(self):
        my_result = self.flux_method.calculate(self.mistis_steps)
        bad_result = {(self.state_A, self.innermost_interface_A): 2.0,
                      (self.state_B, self.innermost_interface_B): 2.0}
        self.flux_method.combine_results(my_result, bad_result)


class TestMinusMoveFlux(TISAnalysisTester):
    def setup(self):
        super(TestMinusMoveFlux, self).setup()

        a = 0.1  # just a number to simplify the trajectory-making
        minus_move_descriptions = [
            [-a, a, a, -a, -a, -a, -a, -a, a, a, a, a, a, -a],
            [-a, a, a, a, -a, -a, -a, a, a, a, -a]
        ]

        engine = RandomMDEngine()  # to get snapshot_timestep

        self.mistis_scheme = paths.DefaultScheme(self.mistis, engine)
        self.mistis_scheme.build_move_decision_tree()
        self.mistis_minus_steps = self._make_fake_minus_steps(
            scheme=self.mistis_scheme,
            descriptions=minus_move_descriptions
        )
        self.mistis_minus_flux = MinusMoveFlux(self.mistis_scheme)

        self.mstis_scheme = paths.DefaultScheme(self.mstis, engine)
        self.mstis_scheme.build_move_decision_tree()
        self.mstis_minus_steps = self._make_fake_minus_steps(
            scheme=self.mstis_scheme,
            descriptions=minus_move_descriptions
        )
        self.mstis_minus_flux = MinusMoveFlux(self.mstis_scheme)

    def _make_fake_minus_steps(self, scheme, descriptions):
        network = scheme.network
        state_adjustment = {
            self.state_A: lambda x: x,
            self.state_B: lambda x: 1.0 - x
        }

        minus_ensemble_to_mover = {m.minus_ensemble: m
                                   for m in scheme.movers['minus']}

        assert_equal(set(minus_ensemble_to_mover.keys()),
                     set(network.minus_ensembles))
        steps = []
        mccycle = 0
        for minus_traj in descriptions:
            for i, minus_ensemble in enumerate(network.minus_ensembles):
                replica = -1 - i
                adjustment = state_adjustment[minus_ensemble.state_vol]
                traj = make_1d_traj([adjustment(s) for s in minus_traj])
                assert_equal(minus_ensemble(traj), True)
                samp = paths.Sample(trajectory=traj,
                                    ensemble=minus_ensemble,
                                    replica=replica)
                sample_set = paths.SampleSet([samp])
                change = paths.AcceptedSampleMoveChange(
                    samples=[samp],
                    mover=minus_ensemble_to_mover[samp.ensemble],
                    details=paths.Details()
                )
                # NOTE: this makes it so that only one ensemble is
                # represented in the same set at any time, which isn't quite
                # how it actually works. However, this is doesn't matter for
                # the current implementation
                steps.append(paths.MCStep(mccycle=mccycle,
                                          active=sample_set,
                                          change=change))

                mccycle += 1
        assert_equal(len(steps), 4)
        return steps

    def test_get_minus_steps(self):
        all_mistis_steps = self.mistis_steps + self.mistis_minus_steps
        mistis_minus_steps = \
            self.mistis_minus_flux._get_minus_steps(all_mistis_steps)
        assert_equal(len(mistis_minus_steps), len(self.mistis_minus_steps))
        assert_items_equal(mistis_minus_steps, self.mistis_minus_steps)
        # this could be repeated for MSTIS, but why?

    def test_calculate(self):
        avg_t_in = (5.0 + 3.0) / 2
        avg_t_out = (2.0 + 5.0 + 3.0 + 3.0) / 4
        expected_flux = 1.0 / (avg_t_in + avg_t_out)

        mistis_flux = \
            self.mistis_minus_flux.calculate(self.mistis_minus_steps)
        for flux in mistis_flux.values():  # all values are the same
            assert_almost_equal(flux, expected_flux)

        mstis_flux = \
            self.mstis_minus_flux.calculate(self.mstis_minus_steps)
        for flux in mstis_flux.values():  # all values are the same
            assert_almost_equal(flux, expected_flux)

    @raises(ValueError)
    def test_bad_network(self):
        # raises error if more than one transition shares a minus ensemble
        # (flux cannot be calculated with multiple interface set minus move)
        state_C = paths.CVDefinedVolume(self.cv_x, 0.5, 0.7)
        trans_AB = self.mistis.transitions[(self.state_A, self.state_B)]
        trans_BA = self.mistis.transitions[(self.state_B, self.state_A)]
        interfaces_AB = trans_AB.interfaces
        interfaces_BA = trans_BA.interfaces
        interfaces_AC = trans_AB.interfaces
        bad_mistis = paths.MISTISNetwork([
            (self.state_A, interfaces_AB, self.state_B),
            (self.state_B, interfaces_BA, self.state_A),
            (self.state_A, interfaces_AC, state_C)
        ])
        scheme = paths.DefaultScheme(bad_mistis)
        scheme.build_move_decision_tree()
        MinusMoveFlux(scheme)


class TestPathLengthHistogrammer(TISAnalysisTester):
    def _check_network_results(self, network, hists):
        results = {0: {(3.0,): 2, (5.0,): 1, (7.0,): 1},
                   1: {(5.0,): 2, (7.0,): 1, (12.0,): 1},
                   2: {(7.0,): 2, (12.0,): 2}}

        ensembles_AB = self.sampling_ensembles_for_transition(network,
                                                              self.state_A,
                                                              self.state_B)
        ensembles_BA = self.sampling_ensembles_for_transition(network,
                                                              self.state_B,
                                                              self.state_A)
        for (key, dct) in results.items():
            hist_dct_AB = hists[ensembles_AB[key]]._histogram
            assert_equal(dict(hist_dct_AB), dct)
            hist_dct_BA = hists[ensembles_BA[key]]._histogram
            assert_equal(dict(hist_dct_BA), dct)

    def test_calculate(self):
        default_histogrammer = \
                PathLengthHistogrammer(self.mistis.sampling_ensembles)
        assert_equal(default_histogrammer.hist_parameters,
                     {'bin_width': 5, 'bin_range': (0, 1000)})

        mistis_histogrammer = PathLengthHistogrammer(
            ensembles=self.mistis.sampling_ensembles,
            hist_parameters={'bin_width': 1, 'bin_range': (0, 10)}
        )
        mistis_hists = mistis_histogrammer.calculate(self.mistis_steps)
        self._check_network_results(self.mistis, mistis_hists)

        mstis_histogrammer = PathLengthHistogrammer(
            ensembles=self.mstis.sampling_ensembles,
            hist_parameters={'bin_width': 1, 'bin_range': (0, 10)}
        )
        mstis_hists = mstis_histogrammer.calculate(self.mstis_steps)
        self._check_network_results(self.mstis, mstis_hists)


class TestFullHistogramMaxLambda(TISAnalysisTester):
    def _check_transition_results(self, transition, hists):
        raw_lambda_results = {
            0: {-0.05: 0, 0.05: 2, 0.15: 1, 0.25: 1},
            1: {-0.05: 0, 0.05: 0, 0.15: 2, 0.25: 1, 0.35: 0, 1.05: 1},
            2: {-0.05: 0, 0.05: 0, 0.15: 0, 0.25: 2, 0.35: 0, 1.05: 2}
        }
        for (ens, result_dct) in raw_lambda_results.items():
            hist = hists[transition.ensembles[ens]]()
            for key in result_dct.keys():
                assert_almost_equal(result_dct[key], hist(key))

    def test_calculate(self):
        mistis_AB = self.mistis.transitions[(self.state_A, self.state_B)]
        mistis_AB_histogrammer = FullHistogramMaxLambdas(
            transition=mistis_AB,
            hist_parameters={'bin_width': 0.1, 'bin_range': (-0.1, 1.1)}
        )
        mistis_AB_hists = mistis_AB_histogrammer.calculate(self.mistis_steps)
        self._check_transition_results(mistis_AB, mistis_AB_hists)

        mistis_BA = self.mistis.transitions[(self.state_B, self.state_A)]
        mistis_BA_histogrammer = FullHistogramMaxLambdas(
            transition=mistis_BA,
            hist_parameters={'bin_width': 0.1, 'bin_range': (-0.1, 1.1)}
        )
        mistis_BA_hists = mistis_BA_histogrammer.calculate(self.mistis_steps)
        self._check_transition_results(mistis_BA, mistis_BA_hists)

        mstis_AB = self.mstis.transitions[(self.state_A, self.state_B)]
        mstis_AB_histogrammer = FullHistogramMaxLambdas(
            transition=mstis_AB,
            hist_parameters={'bin_width': 0.1, 'bin_range': (-0.1, 1.1)}
        )
        mstis_AB_hists = mstis_AB_histogrammer.calculate(self.mstis_steps)
        self._check_transition_results(mstis_AB, mstis_AB_hists)

        mstis_BA = self.mstis.transitions[(self.state_B, self.state_A)]
        mstis_BA_histogrammer = FullHistogramMaxLambdas(
            transition=mstis_BA,
            hist_parameters={'bin_width': 0.1, 'bin_range': (-0.1, 1.1)}
        )
        mstis_BA_hists = mstis_BA_histogrammer.calculate(self.mstis_steps)
        self._check_transition_results(mstis_BA, mstis_BA_hists)

    @raises(RuntimeError)
    def test_calculate_no_max_lambda(self):
        mistis_AB = self.mistis.transitions[(self.state_A, self.state_B)]
        modified_transition = paths.TISTransition(
            stateA=mistis_AB.stateA,
            stateB=mistis_AB.stateB,
            interfaces=mistis_AB.interfaces.volumes,
            orderparameter=mistis_AB.orderparameter
        )
        FullHistogramMaxLambdas(
            transition=modified_transition,
            hist_parameters={'bin_width': 0.1, 'bin_range': (-0.1, 1.1)}
        )


class TestConditionalTransitionProbability(TISAnalysisTester):
    def _check_network_results(self, network, ctp_results):
        results = {0: 0.0, 1: 0.25, 2: 0.5}
        ensembles_AB = self.sampling_ensembles_for_transition(network,
                                                              self.state_A,
                                                              self.state_B)
        ensembles_BA = self.sampling_ensembles_for_transition(network,
                                                              self.state_B,
                                                              self.state_A)

        for ens_num in range(len(ensembles_AB)):
            dct_AB = ctp_results[ensembles_AB[ens_num]]
            result = results[ens_num]
            if result != 0.0:
                assert_equal(dct_AB[self.state_B], result)
            if result != 1.0:
                assert_equal(dct_AB[self.state_A], 1.0-result)

        for ens_num in range(len(ensembles_BA)):
            dct_BA = ctp_results[ensembles_BA[ens_num]]
            result = results[ens_num]
            if result != 0.0:
                assert_equal(dct_BA[self.state_A], result)
            if result != 1.0:
                assert_equal(dct_BA[self.state_B], 1.0-result)

    def test_calculate(self):
        mistis_ctp_calc = ConditionalTransitionProbability(
            ensembles=self.mistis.sampling_ensembles,
            states=[self.state_A, self.state_B]
        )
        mistis_ctp = mistis_ctp_calc.calculate(self.mistis_steps)
        self._check_network_results(self.mistis, mistis_ctp)

        mstis_ctp_calc = ConditionalTransitionProbability(
            ensembles=self.mstis.sampling_ensembles,
            states=[self.state_A, self.state_B]
        )
        mstis_ctp = mstis_ctp_calc.calculate(self.mstis_steps)
        self._check_network_results(self.mstis, mstis_ctp)


class TestTotalCrossingProbability(TISAnalysisTester):
    def test_calculate(self):
        # a bit of integration test, until we make a MaxLambdaStub
        results = {0.0: 1.0, 0.1: 0.5, 0.2: 0.25, 0.3: 0.125,
                   0.5: 0.125, 1.0: 0.125}

        mistis_AB = self.mistis.transitions[(self.state_A, self.state_B)]
        mistis_AB_max_lambda = FullHistogramMaxLambdas(
            transition=mistis_AB,
            hist_parameters={'bin_width': 0.1, 'bin_range': (-0.1, 1.1)}
        )
        mistis_AB_tcp = TotalCrossingProbability(mistis_AB_max_lambda)
        tcp_AB = mistis_AB_tcp.calculate(self.mistis_steps)
        for (x, result) in results.items():
            assert_almost_equal(tcp_AB(x), result)

        mistis_BA = self.mistis.transitions[(self.state_B, self.state_A)]
        mistis_BA_max_lambda = FullHistogramMaxLambdas(
            transition=mistis_BA,
            hist_parameters={'bin_width': 0.1, 'bin_range': (-0.1, 1.1)}
        )
        mistis_BA_tcp = TotalCrossingProbability(mistis_BA_max_lambda)
        tcp_BA = mistis_BA_tcp.calculate(self.mistis_steps)
        for (x, result) in results.items():
            assert_almost_equal(tcp_AB(x), result)

        mstis_AB = self.mstis.transitions[(self.state_A, self.state_B)]
        mstis_AB_max_lambda = FullHistogramMaxLambdas(
            transition=mstis_AB,
            hist_parameters={'bin_width': 0.1, 'bin_range': (-0.1, 1.1)}
        )
        mstis_AB_tcp = TotalCrossingProbability(mstis_AB_max_lambda)
        tcp_AB = mstis_AB_tcp.calculate(self.mstis_steps)
        for (x, result) in results.items():
            assert_almost_equal(tcp_AB(x), result)

        mstis_BA = self.mstis.transitions[(self.state_B, self.state_A)]
        mstis_BA_max_lambda = FullHistogramMaxLambdas(
            transition=mstis_BA,
            hist_parameters={'bin_width': 0.1, 'bin_range': (-0.1, 1.1)}
        )
        mstis_BA_tcp = TotalCrossingProbability(mstis_BA_max_lambda)
        tcp_BA = mstis_BA_tcp.calculate(self.mstis_steps)
        for (x, result) in results.items():
            assert_almost_equal(tcp_BA(x), result)


class TestStandardTransitionProbability(TISAnalysisTester):
    def _check_network_results(self, network, steps):
        for transition in network.transitions.values():
            max_lambda_calc = FullHistogramMaxLambdas(
                transition=transition,
                hist_parameters={'bin_width': 0.1, 'bin_range': (-0.1, 1.1)}
            )
            std_tp = StandardTransitionProbability(
                transition=transition,
                tcp_method=TotalCrossingProbability(max_lambda_calc),
                ctp_method=ConditionalTransitionProbability(
                    ensembles=[transition.ensembles[-1]],
                    states=[self.state_A, self.state_B]
                )
            )
            results = std_tp.calculate(steps)
            assert_almost_equal(results, 0.125)

    def test_calculate(self):
        self._check_network_results(self.mistis, self.mistis_steps)
        self._check_network_results(self.mstis, self.mstis_steps)

    def test_missing_ctp(self):
        ensembles_AB = self.sampling_ensembles_for_transition(
            self.mistis, self.state_A, self.state_B
        )
        ensembles_BA = self.sampling_ensembles_for_transition(
            self.mistis, self.state_B, self.state_A
        )

        all_ensembles = ensembles_AB + ensembles_BA
        replicas = range(len(all_ensembles))
        set_trajectories = [self.trajs_AB[2]]*3 + [self.trajs_BA[2]]*3
        zipped = list(zip(set_trajectories, all_ensembles, replicas))
        mover_stub_mistis = MoverWithSignature(self.mistis.all_ensembles,
                                               self.mistis.all_ensembles)
        sample_sets = []
        n_steps = 3
        for step in range(n_steps):
            sample_set = paths.SampleSet([
                paths.Sample(trajectory=traj,
                             ensemble=ens,
                             replica=rep)
                for (traj, ens, rep) in zipped
            ])
            sample_set.sanity_check()
            sample_sets.append(sample_set)

        steps = self._make_fake_steps(sample_sets, mover_stub_mistis)
        for transition in self.mistis.transitions.values():
            max_lambda_calc = FullHistogramMaxLambdas(
                transition=transition,
                hist_parameters={'bin_width': 0.1, 'bin_range': (-0.1, 1.1)}
            )
            std_tp = StandardTransitionProbability(
                transition=transition,
                tcp_method=TotalCrossingProbability(max_lambda_calc),
                ctp_method=ConditionalTransitionProbability(
                    ensembles=[transition.ensembles[-1]],
                    states=[self.state_A, self.state_B]
                )
            )
            results = std_tp.calculate(steps)
            assert_almost_equal(results, 0.0)


class TestTransitionDictResults(TISAnalysisTester):
    def setup(self):
        super(TestTransitionDictResults, self).setup()
        results_dict = {(self.state_A, self.state_B): 1,
                        (self.state_B, self.state_A): 2}
        self.mistis_transition_dict = TransitionDictResults(
            results_dict=results_dict,
            network=self.mistis,
            allow_sampling=False
        )
        self.mstis_transition_dict = TransitionDictResults(
            results_dict=results_dict,
            network=self.mstis,
            allow_sampling=True
        )

    def test_iter(self):
        assert_equal(set(pair for pair in self.mistis_transition_dict),
                     set(pair for pair in self.mstis_transition_dict))

    def test_get_by_pair(self):
        assert_equal(
            self.mstis_transition_dict[(self.state_A, self.state_B)], 1
        )
        assert_equal(
            self.mstis_transition_dict[(self.state_A, self.state_B)],
            self.mistis_transition_dict[(self.state_A, self.state_B)]
        )
        assert_equal(
            self.mstis_transition_dict[(self.state_B, self.state_A)], 2
        )
        assert_equal(
            self.mstis_transition_dict[(self.state_B, self.state_A)],
            self.mistis_transition_dict[(self.state_B, self.state_A)]
        )

    @raises(KeyError)
    def test_get_bad_pair(self):
        self.mistis_transition_dict[(self.state_A, self.state_A)]

    def test_get_by_transition(self):
        mistis_AB = self.mistis.transitions[(self.state_A, self.state_B)]
        mstis_AB = self.mstis.transitions[(self.state_A, self.state_B)]
        assert_equal(self.mistis_transition_dict[mistis_AB], 1)
        assert_equal(self.mistis_transition_dict[mistis_AB],
                     self.mstis_transition_dict[mstis_AB])

    def test_get_by_sampling_transition(self):
        from_A = self.mstis.from_state[self.state_A]
        from_B = self.mstis.from_state[self.state_B]
        assert_equal(self.mstis_transition_dict[from_A], 1)
        assert_equal(self.mstis_transition_dict[from_B], 2)

    @raises(KeyError)
    def test_bad_get_sampling_transition(self):
        sampling_trans = self.mistis.sampling_transitions[0]
        self.mistis_transition_dict[sampling_trans]

    def test_to_pandas(self):
        result = [[float("nan"), 1], [2, float("nan")]]
        order = [self.state_A, self.state_B]
        ordered_names = ["A", "B"]
        pd_result = pd.DataFrame(data=result, index=ordered_names,
                                 columns=ordered_names)
        pd_mistis = self.mistis_transition_dict.to_pandas(order=order)
        pd_mstis = self.mstis_transition_dict.to_pandas()
        assert_frame_equal(pd_mistis, pd_mstis)
        assert_frame_equal(pd_mistis, pd_result)


class TestTISAnalysis(TISAnalysisTester):
    def _make_tis_analysis(self, network):
        # NOTE: this might be useful as a description of the overall nested
        # structure of these TISAnalysis objects. However, this isn't a
        # practical way to do the calculation, because it repeats effort
        # (hence the StandardTISAnalysis object)
        tis_analysis = TISAnalysis(
            network=network,
            flux_method=DictFlux({
                (t[0].stateA, t[0].interfaces[0]): 0.1
                for t in network.special_ensembles['minus'].values()
            }),
            transition_probability_methods={
                transition: StandardTransitionProbability(
                    transition=transition,
                    tcp_method=TotalCrossingProbability(
                        max_lambda_calc=FullHistogramMaxLambdas(
                            transition=transition,
                            hist_parameters={'bin_width': 0.1,
                                             'bin_range': (-0.1, 1.1)}
                        )
                    ),
                    ctp_method=ConditionalTransitionProbability(
                        ensembles=[transition.ensembles[-1]],
                        states=[self.state_A, self.state_B]
                    )
                )
                for transition in network.transitions.values()
            }
        )
        return tis_analysis

    def setup(self):
        super(TestTISAnalysis, self).setup()
        self.mistis_analysis = self._make_tis_analysis(self.mistis)
        self.mistis_analysis.calculate(self.mistis_steps)
        self.mstis_analysis = self._make_tis_analysis(self.mstis)
        self.mstis_analysis.calculate(self.mstis_steps)

    def test_bad_access_cached_results(self):
        no_results = self._make_tis_analysis(self.mistis)
        _ = self.mistis_analysis._access_cached_result('rate')
        # use a try/except here instead of @raises so that we also test that
        # the calculated version (previous line) works as expected
        try:
            no_results._access_cached_result('rate')
        except AttributeError:
            pass  # this is the expected test result

    def test_flux_matrix(self):
        assert_equal(self.mistis_analysis.flux_matrix,
                     {(t.stateA, t.interfaces[0]): 0.1
                      for t in self.mistis.sampling_transitions})
        assert_equal(self.mstis_analysis.flux_matrix,
                     {(t.stateA, t.interfaces[0]): 0.1
                      for t in self.mstis.sampling_transitions})

    def test_flux(self):
        for transition in self.mistis.sampling_transitions:
            state = transition.stateA
            innermost = transition.interfaces[0]
            assert_equal(self.mistis_analysis.flux(state, innermost), 0.1)

        for transition in self.mstis.sampling_transitions:
            state = transition.stateA
            innermost = transition.interfaces[0]
            assert_equal(self.mstis_analysis.flux(state, innermost), 0.1)

    def test_flux_through_state(self):
        flux_dict = {(t.stateA, t.interfaces[0]): 0.1
                     for t in self.mistis.sampling_transitions}
        flux_dict.update({(self.state_A, self.state_A): 0.5})
        tis = TISAnalysis(
            network=self.mistis,
            flux_method=DictFlux(flux_dict),
            transition_probability_methods={
                transition: StandardTransitionProbability(
                    transition=transition,
                    tcp_method=TotalCrossingProbability(
                        max_lambda_calc=FullHistogramMaxLambdas(
                            transition=transition,
                            hist_parameters={'bin_width': 0.1,
                                             'bin_range': (-0.1, 1.1)}
                        )
                    ),
                    ctp_method=ConditionalTransitionProbability(
                        ensembles=[transition.ensembles[-1]],
                        states=[self.state_A, self.state_B]
                    )
                )
                for transition in self.mistis.transitions.values()
            }
        )
        tis.calculate(self.mistis_steps)
        trans_AB = self.mistis.transitions[(self.state_A, self.state_B)]
        assert_equal(tis.flux(self.state_A, trans_AB.interfaces[0]), 0.1)
        assert_equal(tis.flux(self.state_A), 0.5)

    def test_state_fluxes(self):
        for transition in self.mistis.sampling_transitions:
            state = transition.stateA
            innermost = transition.interfaces[0]
            assert_equal(self.mistis_analysis.state_fluxes(state),
                         {(state, innermost): 0.1})

        for transition in self.mstis.sampling_transitions:
            state = transition.stateA
            innermost = transition.interfaces[0]
            assert_equal(self.mstis_analysis.state_fluxes(state),
                         {(state, innermost): 0.1})

    def test_transition_probability_matrix(self):
        pairs = [(self.state_A, self.state_B), (self.state_B, self.state_A)]
        mistis_tp = self.mistis_analysis.transition_probability_matrix
        mstis_tp = self.mstis_analysis.transition_probability_matrix
        for trans_pair in pairs:
            assert_almost_equal(mistis_tp[trans_pair], 0.125)
            assert_almost_equal(mstis_tp[trans_pair], 0.125)

    def test_transition_probability(self):
        pairs = [(self.state_A, self.state_B), (self.state_B, self.state_A)]
        for (vol_1, vol_2) in pairs:
            assert_almost_equal(
                self.mistis_analysis.transition_probability(vol_1, vol_2),
                0.125
            )
            assert_almost_equal(
                self.mstis_analysis.transition_probability(vol_1, vol_2),
                0.125
            )

    def test_rate_matrix(self):
        pairs = [(self.state_A, self.state_B), (self.state_B, self.state_A)]
        mistis_rate = self.mistis_analysis.rate_matrix()
        mstis_rate = self.mstis_analysis.rate_matrix()
        for (vol_1, vol_2) in pairs:
            assert_almost_equal(mistis_rate[(vol_1, vol_2)], 0.0125)
            assert_almost_equal(mstis_rate[(vol_1, vol_2)], 0.0125)

    def test_rate_matrix_calculation(self):
        mistis_analysis = self._make_tis_analysis(self.mistis)
        mistis_rate = mistis_analysis.rate_matrix(steps=self.mistis_steps)
        pairs = [(self.state_A, self.state_B), (self.state_B, self.state_A)]
        for (vol_1, vol_2) in pairs:
            assert_almost_equal(mistis_rate[(vol_1, vol_2)], 0.0125)

    def test_rate(self):
        pairs = [(self.state_A, self.state_B), (self.state_B, self.state_A)]
        for (vol_1, vol_2) in pairs:
            assert_almost_equal(self.mistis_analysis.rate(vol_1, vol_2),
                                0.0125)
            assert_almost_equal(self.mstis_analysis.rate(vol_1, vol_2),
                                0.0125)


class TestStandardTISAnalysis(TestTISAnalysis):
    # inherit from TestTISAnalysis to retest all the same results
    def _make_tis_analysis(self, network, steps=None):
        tis_analysis = StandardTISAnalysis(
            network=network,
            flux_method=DictFlux({(t.stateA, t.interfaces[0]): 0.1
                                  for t in network.sampling_transitions}),
            max_lambda_calcs={t: {'bin_width': 0.1,
                                  'bin_range': (-0.1, 1.1)}
                              for t in network.sampling_transitions},
            steps=steps
        )
        return tis_analysis

    def test_init_with_steps(self):
        mistis_analysis = self._make_tis_analysis(self.mistis,
                                                  self.mistis_steps)
        pairs = [(self.state_A, self.state_B), (self.state_B, self.state_A)]
        for (vol_1, vol_2) in pairs:
            assert_almost_equal(
                mistis_analysis.rate_matrix()[(vol_1, vol_2)],
                0.0125
            )

    def test_crossing_probability(self):
        results = {
            0: {0.0: 1.0, 0.1: 0.5, 0.2: 0.25},
            1: {0.0: 1.0, 0.1: 1.0, 0.2: 0.5, 0.3: 0.25, 1.0: 0.25},
            2: {0.0: 1.0, 0.1: 1.0, 0.2: 1.0, 0.3: 0.5, 1.0: 0.5}
        }

        def check_cp(transition, analysis, results):
            # a little nested function to that does the actual check
            for (i, ens) in enumerate(transition.ensembles):
                cp_ens = analysis.crossing_probability(ens)
                for x in results[i]:
                    assert_almost_equal(results[i][x], cp_ens(x))

        for (network, analysis) in [(self.mistis, self.mistis_analysis),
                                    (self.mstis, self.mstis_analysis)]:
            for transition in network.transitions.values():
                check_cp(transition, analysis, results)

    def test_conditional_transition_probability(self):
        expected_data = [[0.5, 0.5], [0.5, 0.5]]
        states = ['A', 'B']
        mistis_interfaces = ['A->B 2', 'B->A 2']
        mstis_interfaces = ['Out A 2', 'Out B 2']

        expected_mistis = pd.DataFrame(data=expected_data,
                                       index=mistis_interfaces,
                                       columns=states)
        mistis_ctp = self.mistis_analysis.conditional_transition_probability
        assert_equal(set(states), set(mistis_ctp.columns))
        assert_equal(set(mistis_interfaces), set(mistis_ctp.index))
        for iface in mistis_interfaces:
            for state in states:
                assert_equal(expected_mistis.loc[(iface, state)],
                             mistis_ctp.loc[(iface, state)])

        expected_mstis = pd.DataFrame(data=expected_data,
                                      index=mstis_interfaces,
                                      columns=states)
        mstis_ctp = self.mstis_analysis.conditional_transition_probability
        assert_equal(set(states), set(mstis_ctp.columns))
        assert_equal(set(mstis_interfaces), set(mstis_ctp.index))
        for iface in mstis_interfaces:
            for state in states:
                assert_equal(expected_mstis.loc[(iface, state)],
                             mstis_ctp.loc[(iface, state)])

    def test_total_crossing_probability(self):
        results = {0.0: 1.0, 0.1: 0.5, 0.2: 0.25, 0.3: 0.125,
                   0.5: 0.125, 1.0: 0.125}

        mistis_tcp = self.mistis_analysis.total_crossing_probability
        for transition in self.mistis.transitions.values():
            tcp = mistis_tcp[transition]
            for x in results:
                assert_almost_equal(results[x], tcp(x))

        mstis_tcp = self.mstis_analysis.total_crossing_probability
        for transition in self.mstis.transitions.values():
            tcp = mstis_tcp[transition]
            for x in results:
                assert_almost_equal(results[x], tcp(x))

    @raises(TypeError)
    def test_bad_no_flux(self):
        network = self.mistis
        StandardTISAnalysis(
            network=network,
            max_lambda_calcs={t: {'bin_width': 0.1,
                                  'bin_range': (-0.1, 1.1)}
                              for t in network.sampling_transitions}
        )

    @raises(RuntimeError)
    def test_bad_max_lambda_calcs(self):
        network = self.mistis
        StandardTISAnalysis(
            network=network,
            flux_method=DictFlux({(t.stateA, t.interfaces[0]): 0.1
                                  for t in network.sampling_transitions})
        )

    def test_init_ensemble_histogrammer_max_lambda(self):
        network = self.mistis
        max_lambda_calcs = {
            t: FullHistogramMaxLambdas(
                transition=t,
                hist_parameters={'bin_width': 0.1, 'bin_range': (-0.1, 1.1)}
            )
            for t in network.sampling_transitions
        }
        tis_analysis = StandardTISAnalysis(
            network=network,
            flux_method=DictFlux({(t.stateA, t.interfaces[0]): 0.1
                                  for t in network.sampling_transitions}),
            max_lambda_calcs=max_lambda_calcs,
            steps=self.mistis_steps
        )
        rate = tis_analysis.rate_matrix()
        pairs = [(self.state_A, self.state_B), (self.state_B, self.state_A)]
        for (vol_1, vol_2) in pairs:
            assert_almost_equal(rate[(vol_1, vol_2)], 0.0125)

    def test_with_minus_move_flux(self):
        network = self.mstis
        scheme = paths.DefaultScheme(network, engine=RandomMDEngine())
        scheme.build_move_decision_tree()

        # create the minus move steps
        # `center` is the edge of the state/innermost interface
        center = {self.state_A: 0.0, self.state_B: 1.0}
        replica = {self.state_A: -1, self.state_B: -2}
        minus_ensemble_to_mover = {m.minus_ensemble: m
                                   for m in scheme.movers['minus']}
        state_to_minus_ensemble = {ens.state_vol: ens
                                   for ens in network.minus_ensembles}
        minus_changes = []
        # `delta` is the change on either side for in vs. out
        for (state, delta) in [(self.state_A, 0.1), (self.state_B, -0.1)]:
            minus_ens = state_to_minus_ensemble[state]
            minus_mover = minus_ensemble_to_mover[minus_ens]
            a_in = center[state] - delta
            a_out = center[state] + delta
            # note that these trajs are equivalent to minus move
            # descriptions in TestMinusMoveFlux
            seq_1 = [a_in] + [a_out]*2 + [a_in]*5 + [a_out]*5 + [a_in]
            seq_2 = [a_in] + [a_out]*3 + [a_in]*3 + [a_out]*3 + [a_in]

            for seq in [seq_1, seq_2]:
                traj = make_1d_traj(seq)
                assert_equal(minus_ens(traj), True)
                samp = paths.Sample(trajectory=traj,
                                    ensemble=minus_ens,
                                    replica=replica[state])
                _ = paths.SampleSet([samp])
                change = paths.AcceptedSampleMoveChange(
                    samples=[samp],
                    mover=minus_mover,
                    details=paths.Details()
                )
                minus_changes.append(change)

        active = self.mstis_steps[0].active
        steps = []
        cycle = -1
        for m_change in minus_changes:
            cycle += 1
            active = active.apply_samples(m_change.samples)
            step = paths.MCStep(mccycle=cycle,
                                active=active,
                                change=m_change)
            steps.append(step)
            for old_step in self.mstis_steps[1:]:
                cycle += 1
                active = active.apply_samples(old_step.change.samples)
                step = paths.MCStep(mccycle=cycle,
                                    active=active,
                                    change=old_step.change)
                steps.append(step)

        analysis = StandardTISAnalysis(
            network=self.mstis,
            scheme=scheme,
            max_lambda_calcs={t: {'bin_width': 0.1,
                                  'bin_range': (-0.1, 1.1)}
                              for t in network.sampling_transitions},
            steps=steps
        )

        # now we actually verify correctness
        avg_t_in = (5.0 + 3.0) / 2
        avg_t_out = (2.0 + 5.0 + 3.0 + 3.0) / 4
        expected_flux = 1.0 / (avg_t_in + avg_t_out)

        # NOTE: Apparently this approach screws up the TCP calculation. I
        # think this is a problem in the fake data, not the simulation.
        for flux in analysis.flux_matrix.values():
            assert_almost_equal(flux, expected_flux)

    @pytest.mark.parametrize('progress', ['all', 'default', 'none',
                                          'tqdm', 'silent'])
    def test_progress_setter(self, progress):
        analysis = self.mstis_analysis
        analysis.progress = progress
        expected_flux, expected_ctp, expected_max_lambda = {
            'all': (True, True, True),
            'default': (True, True, False),
            'none': (False, False, False),
            'tqdm': (True, True, False),
            'silent': (True, True, False),
        }[progress]
        flux_method = analysis.flux_method
        assert flux_method.progress.keywords['leave'] is expected_flux
        ctp_method = analysis.ctp_method
        assert ctp_method.progress.keywords['leave'] is expected_ctp
        max_lambda_methods = [tcp.max_lambda_calc
                              for tcp in analysis.tcp_methods.values()]
        for max_lambda in max_lambda_methods:
            prog = max_lambda.progress
            assert prog.keywords['leave'] is expected_max_lambda




