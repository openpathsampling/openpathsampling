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
        transition = self.network.transitions[(self.stateA, self.stateB)]
        ensembles = transition.ensembles
        # create the biases
        bias_table = {}
        bias_table[ensembles[0]] = 1.0
        bias_table[ensembles[1]] = 0.5
        bias_table[ensembles[2]] = 0.2
        self.bias = BiasEnsembleTable(bias_table)
        # samples, moves, changes
        traj = make_1d_traj([-0.55, -0.45, -0.35, -0.25, -0.15, -0.26,
                             -0.36, -0.46, -0.56])
        s0 = paths.Sample(replica=0, ensemble=ensembles[0], trajectory=traj)
        s1 = paths.Sample(replica=1, ensemble=ensembles[1], trajectory=traj)
        s2 = paths.Sample(replica=2, ensemble=ensembles[2], trajectory=traj)
        self.sample_set = paths.SampleSet([s0, s1, s2])
        self.sample_set.sanity_check()
        move_01 = paths.EnsembleHopMover(ensembles[0], ensembles[1])
        move_02 = paths.EnsembleHopMover(ensembles[0], ensembles[2])
        move_12 = paths.EnsembleHopMover(ensembles[1], ensembles[2])
        move_21 = paths.EnsembleHopMover(ensembles[2], ensembles[1])
        move_20 = paths.EnsembleHopMover(ensembles[2], ensembles[0])
        move_10 = paths.EnsembleHopMover(ensembles[1], ensembles[0])
        # NOTE: all changes here are accepted
        self.change_01 = move_01.move(self.sample_set)
        self.change_02 = move_02.move(self.sample_set)
        self.change_12 = move_12.move(self.sample_set)
        self.change_21 = move_21.move(self.sample_set)
        self.change_20 = move_20.move(self.sample_set)
        self.change_10 = move_10.move(self.sample_set)
        # convenience lists for changes going outward vs. inward
        self.out_changes = [self.change_01, self.change_02, self.change_12]
        self.in_changes = [self.change_10, self.change_20, self.change_21]


    def test_bias_ensemble_new_to_old(self):
        # The o->n change is the denominator of the ratio.

        # for old_to_new, the probability of moving outerward depends on the
        # ratio of the probabilities of the two ensembles
        change_vals = { 
            self.change_01 : 0.5,
            self.change_02 : 0.2,
            self.change_12 : 0.2 / 0.5,
            self.change_10 : 1.0 / 0.5,
            self.change_20 : 1.0 / 0.2,
            self.change_21 : 0.5 / 0.2
        }
        for change in change_vals.keys():
            test_val = min(1.0, change_vals[change])
            assert_almost_equal(
                self.bias.probability_new_to_old(self.sample_set, change),
                test_val
            )

    def test_bias_ensemble_old_to_new(self):
        # The n->o change is the numerator of the ratio.

        # prob of moving inward is the ratio of the interface weights (cf
        # test_bias_ensemble_old_to_new)
        change_vals = {
            self.change_10 : 0.5,
            self.change_20 : 0.2,
            self.change_21 : 0.2 / 0.5,
            self.change_01 : 1.0 / 0.5,
            self.change_02 : 1.0 / 0.2,
            self.change_12 : 0.5 / 0.2
        }
        for change in change_vals.keys():
            test_val = min(1.0, change_vals[change])
            assert_almost_equal(
                self.bias.probability_old_to_new(self.sample_set, change),
                test_val
            )

    def test_combo_bias(self):
        #TODO: we don't support this yet, but we need to at some point
        # test what happens if you have more than one sample in the change
        # if you do a multi-step bias (in the same direction), then it seems
        # to me that the total bias should be the same as if it were a
        # one-step. This is *not* true if hops are in different directions.
        # Then it depends on the product of the "upstream" hops

        ensembles = self.network.transitions[(self.stateA, self.stateB)].ensembles
        # all downstream move
        move_012 = paths.SequentialMover([
            paths.EnsembleHopMover(ensembles[0], ensembles[1]),
            paths.EnsembleHopMover(ensembles[1], ensembles[2])
        ])
        change_012 = move_012.move(self.sample_set)

        # all upstream move
        move_210 = paths.SequentialMover([
            paths.EnsembleHopMover(ensembles[2], ensembles[1]),
            paths.EnsembleHopMover(ensembles[1], ensembles[0])
        ])
        change_210 = move_210.move(self.sample_set)

        # assert_almost_equal(
            # self.bias.probability_old_to_new(change_210, self.sample_set), 1.0
        # )
        # assert_almost_equal(
            # self.bias.probability_new_to_old(change_210, self.sample_set), 0.2
        # )
        raise SkipTest


