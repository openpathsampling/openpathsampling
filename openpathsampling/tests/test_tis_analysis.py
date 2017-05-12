from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        assert_almost_equal, raises)
from nose.plugins.skip import Skip, SkipTest
from test_helpers import (
    true_func, assert_equal_array_array, make_1d_traj, data_filename,
    MoverWithSignature
)

from openpathsampling.analysis.tis_analysis import *

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
        analysis_AB = network.transitions[(self.state_A, self.state_B)]
        analysis_BA = network.transitions[(self.state_B, self.state_A)]

        sampling_AB = network.analysis_to_sampling[analysis_AB][0]
        sampling_BA = network.analysis_to_sampling[analysis_BA][0]

        all_ensembles = sampling_AB.ensembles + sampling_BA.ensembles
        replicas = range(len(all_ensembles))

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
            zipped = zip(set_trajectories, all_ensembles, replicas)
            sample_set = paths.SampleSet([
                paths.Sample(trajectory=traj,
                             ensemble=ens,
                             replica=rep)
                for (traj, ens, rep) in zip(set_trajectories, all_ensembles,
                                            range(len(all_ensembles)))
            ])
            sample_sets.append(sample_set)
        return sample_sets

    def setup(self):
        # set up the trajectories, ensembles, etc. for this test
        cv_A = paths.FunctionCV('Id', lambda s: s.xyz[0][0])
        cv_B = paths.FunctionCV('-Id', lambda s: 1.0-s.xyz[0][0])
        self.state_A = paths.CVDefinedVolume(cv_A, float("-inf"), 0.0)
        self.state_B = paths.CVDefinedVolume(cv_B, float("-inf"), 0.0)
        interfaces_AB = paths.VolumeInterfaceSet(cv_A, float("-inf"),
                                                 [0.0, 0.1, 0.2])
        interfaces_BA = paths.VolumeInterfaceSet(cv_B, float("-inf"),
                                                 [0.0, 0.1, 0.2])

        # trajectory that crosses each interface, one state-to-state
        self.trajs_AB = [make_tis_traj_fixed_steps(i) for i in [0, 1, 2]]
        self.trajs_AB += make_1d_traj([(-0.5 + i) * 0.1 for i in range(12)])

        self.trajs_BA = [make_tis_traj_fixed_steps(i, reverse=True)
                         for i in [0, 1, 2]]
        self.trajs_BA += make_1d_traj([1.0 - (-0.5 + i) * 0.1
                                       for i in range(12)])

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
        pass


class TestWeightedTrajectories(TISAnalysisTester):
    def test_steps_to_weighted_trajectories(self):
        raise SkipTest


class TestDictFlux(TISAnalysisTester):
    def setup(self):
        super(TestDictFlux, self).setup()
        self.flux_dict = {}  # TODO: add actual values here
        self.flux_method = DictFlux(self.flux_dict)

    def test_calculate(self):
        assert_equal(self.flux_method.calculate(self.mistis_steps),
                     self.flux_dict)

    def test_from_weighted_trajectories(self):
        raise SkipTest


class TestPathLengthHistogrammer(TISAnalysisTester):
    def test_calculate(self):
        raise SkipTest

    def test_from_weighted_trajectories(self):
        raise SkipTest


class TestFullHistogramMaxLambda(TISAnalysisTester):
    def test_calculate(self):
        raise SkipTest

    def test_from_weighted_trajectories(self):
        raise SkipTest


class TestTotalCrossingProbability(TISAnalysisTester):
    def test_calculate(self):
        raise SkipTest

    def test_from_weighted_trajectories(self):
        raise SkipTest
