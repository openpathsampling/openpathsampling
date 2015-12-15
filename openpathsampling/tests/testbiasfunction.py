from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        assert_almost_equal, raises)

from nose.plugins.skip import Skip, SkipTest
from test_helpers import true_func, assert_equal_array_array, make_1d_traj

import openpathsampling as paths
from openpathsampling.bias_function import *

import logging
from openpathsampling import VolumeFactory as vf

quiet_loggers = ["initialization", "ensemble", "netcdfplus.objects",
                 "netcdfplus.netcdfplus", "pathmover"]
for logger in quiet_loggers:
    logging.getLogger("openpathsampling."+logger).setLevel(logging.CRITICAL)


class testBiasEnsembleTable(object):
    def setup(self):
        # create the network
        xval = paths.CV_Function(name="xA", f=lambda s : s.xyz[0][0])
        self.stateA = paths.CVRangeVolume(xval, float("-inf"), -0.5)
        self.stateB = paths.CVRangeVolume(xval, 0.5, float("inf"))
        ifacesA = vf.CVRangeVolumeSet(xval, float("-inf"), 
                                      [-0.5, -0.4, -0.3, -0.2])
        self.network = paths.MISTISNetwork([
            (self.stateA, ifacesA, xval, self.stateB)
        ])
        # create the biases
        transition = self.network.transitions[(self.stateA, self.stateB)]
        ensembles = transition.ensembles
        bias_table = {}
        bias_table[ensembles[0]] = 1.0
        bias_table[ensembles[1]] = 0.5
        bias_table[ensembles[2]] = 0.2
        self.bias = BiasEnsembleTable(bias_table)
        traj = make_1d_traj([-0.55, -0.45, -0.35, -0.25, -0.15, -0.26,
                             -0.36, -0.46, -0.56])
        s0 = paths.Sample(replica=0, ensemble=ensembles[0], trajectory=traj)
        s1 = paths.Sample(replica=1, ensemble=ensembles[1], trajectory=traj)
        s2 = paths.Sample(replica=2, ensemble=ensembles[2], trajectory=traj)
        self.sample_set = paths.SampleSet([s0, s1, s2])
        move_01 = paths.EnsembleHopMover(ensembles[0], ensembles[1])
        move_02 = paths.EnsembleHopMover(ensembles[0], ensembles[2])
        move_12 = paths.EnsembleHopMover(ensembles[1], ensembles[2])
        move_21 = paths.EnsembleHopMover(ensembles[2], ensembles[1])
        move_20 = paths.EnsembleHopMover(ensembles[2], ensembles[0])
        move_10 = paths.EnsembleHopMover(ensembles[1], ensembles[0])
        self.change_01 = move_01.move(self.sample_set)
        self.change_02 = move_02.move(self.sample_set)
        self.change_12 = move_12.move(self.sample_set)
        self.change_21 = move_21.move(self.sample_set)
        self.change_20 = move_20.move(self.sample_set)
        self.change_10 = move_10.move(self.sample_set)
        # convenience lists for changes going outward vs. inward
        self.out_changes = [self.change_01, self.change_02, self.change_12]
        self.in_changes = [self.change_10, self.change_20, self.change_21]


    def test_bias_ensemble_old_to_new(self):
        # The o->n change is the denominator of the ratio.
        for change in self.in_changes:
            assert_almost_equal(
                self.bias.probability_old_to_new(self.sample_set, change),
                1.0
            )

        change_vals = { 
            self.change_01 : 0.5,
            self.change_02 : 0.2,
            self.change_12 : 0.2 / 0.5
        }
        for change in change_vals.keys():
            assert_almost_equal(
                self.bias.probability_old_to_new(self.sample_set, change),
                change_vals[change]
            )
        pass

    def test_bias_ensemble_new_to_old(self):
        pass

    def test_bias_ensemble_prob_ratio(self):
        pass


