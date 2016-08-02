from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        assert_almost_equal, raises, assert_in)

from nose.plugins.skip import Skip, SkipTest
from test_helpers import true_func, assert_equal_array_array, make_1d_traj

import numpy as np

import openpathsampling as paths
from openpathsampling.bias_function import *

import logging
from openpathsampling import VolumeFactory as vf

quiet_loggers = ["initialization", "ensemble", "netcdfplus.objects",
                 "netcdfplus.netcdfplus", "pathmover", "netcdfplus.base"]
for logger in quiet_loggers:
    logging.getLogger("openpathsampling."+logger).setLevel(logging.CRITICAL)


class testBiasEnsembleTable(object):
    def setup(self):
        # create the network
        xval = paths.CV_Function(name="xA", f=lambda s : s.xyz[0][0])
        self.stateA = paths.CVRangeVolume(xval, -1.0, -0.5).named("A")
        self.stateB = paths.CVRangeVolume(xval, 0.5, float("inf")).named("B")
        ifacesA = paths.VolumeInterfaceSet(xval, float(-1.0), 
                                           [-0.5, -0.4, -0.3, -0.2])
        self.network = paths.MISTISNetwork([
            (self.stateA, ifacesA, self.stateB)
        ])
        transition = self.network.transitions[(self.stateA, self.stateB)]
        ensembles = transition.ensembles
        self.xval = xval
        self.ifacesA = ifacesA
        # create the biases
        bias_table = {}
        bias_table[ensembles[0]] = 1.0
        bias_table[ensembles[1]] = 0.5
        bias_table[ensembles[2]] = 0.2
        self.bias = BiasEnsembleTable.ratios_from_dictionary(bias_table)
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

    def test_add_biases(self):
        # this is where we combine multiple biases into one
        xval2 = paths.CV_Function(name="xB", f=lambda s : 0.5-s.xyz[0][0])
        ifacesB = paths.VolumeInterfaceSet(xval2, float("-inf"),
                                           [0.0, 0.1, 0.2, 0.3])
        xval3 = paths.CV_Function(name="xC", f=lambda s : s.xyz[0][0]-2.0)
        stateC = paths.CVRangeVolume(self.xval, -3.0, 2.0)
        ifacesC = paths.VolumeInterfaceSet(xval3, -1.0, [0.0, 0.1, 0.2, 0.3])
        network = paths.MISTISNetwork([
            (self.stateA, self.ifacesA, self.stateB),
            (self.stateB, ifacesB, self.stateA),
            (stateC, ifacesC, self.stateA)
        ])
        ens_A = network.transitions[(self.stateA, self.stateB)].ensembles
        ens_B = network.transitions[(self.stateB, self.stateA)].ensembles
        ens_C = network.transitions[(stateC, self.stateA)].ensembles
        ms_outer = network.special_ensembles['ms_outer'].keys()[0]
        dict_A = {ens_A[0]: 1.0,
                  ens_A[1]: 0.5,
                  ens_A[2]: 0.2,
                  ms_outer: 0.1}
        dict_B = {ens_B[0]: 1.0,
                  ens_B[1]: 0.6,
                  ens_B[2]: 0.3,
                  ms_outer: 0.15}
        dict_C = {ens_C[0]: 1.0,
                  ens_C[1]: 0.8,
                  ens_C[2]: 0.2}

        bias_A = BiasEnsembleTable.ratios_from_dictionary(dict_A)
        bias_B = BiasEnsembleTable.ratios_from_dictionary(dict_B)
        bias_C = BiasEnsembleTable.ratios_from_dictionary(dict_C)
        bias_AB = bias_A + bias_B
        # check the ensembles_to_ids
        assert_equal(len(bias_AB.ensembles_to_ids), 7)
        for ens in ens_A:
            assert_in(bias_AB.ensembles_to_ids[ens], [0, 1, 2])
        for ens in ens_B:
            assert_in(bias_AB.ensembles_to_ids[ens], [3, 4, 5])
        assert_equal(bias_AB.ensembles_to_ids[ms_outer], 6)

        # check values
        df_A = bias_A.dataframe
        df_B = bias_B.dataframe
        df_AB = bias_AB.dataframe
        col_A_msouter = bias_A.ensembles_to_ids[ms_outer]
        col_B_msouter = bias_B.ensembles_to_ids[ms_outer]
        col_AB_msouter = bias_AB.ensembles_to_ids[ms_outer]

        for ens1 in ens_A:
            idx_A = bias_A.ensembles_to_ids[ens1]
            idx_AB = bias_AB.ensembles_to_ids[ens1]
            for ens2 in ens_A:
                col_A = bias_A.ensembles_to_ids[ens2]
                col_AB = bias_AB.ensembles_to_ids[ens2]
                val_A = df_A.loc[idx_A, col_A]
                val_AB = df_AB.loc[idx_AB, col_AB]
                assert_equal(val_A, val_AB)
            for ens2 in ens_B:
                col_AB = bias_AB.ensembles_to_ids[ens2]
                assert_equal(np.isnan(df_AB.loc[idx_AB, col_AB]), True)
            assert_equal(df_A.loc[idx_A, col_A_msouter],
                         df_AB.loc[idx_AB, col_AB_msouter])
            assert_equal(df_A.loc[col_A_msouter, idx_A],
                         df_AB.loc[col_AB_msouter, idx_AB])

        for ens1 in ens_B:
            idx_B = bias_B.ensembles_to_ids[ens1]
            idx_AB = bias_AB.ensembles_to_ids[ens1]
            for ens2 in ens_B:
                col_B = bias_B.ensembles_to_ids[ens2]
                col_AB = bias_AB.ensembles_to_ids[ens2]
                val_B = df_B.loc[idx_B, col_B]
                val_AB = df_AB.loc[idx_AB, col_AB]
                assert_equal(val_B, val_AB)
            for ens2 in ens_A:
                col_AB = bias_AB.ensembles_to_ids[ens2]
                assert_equal(np.isnan(df_AB.loc[idx_AB, col_AB]), True)
            assert_equal(df_B.loc[idx_B, col_B_msouter],
                         df_AB.loc[idx_AB, col_AB_msouter])
            assert_equal(df_B.loc[col_B_msouter, idx_B],
                         df_AB.loc[col_AB_msouter, idx_AB])

        # just to make sure no errors raise when there are NaNs in table
        bias_ABC = bias_A + bias_B + bias_C

