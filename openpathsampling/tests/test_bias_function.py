from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from builtins import range
from builtins import object
from past.utils import old_div
from nose.tools import (assert_equal, assert_not_equal, assert_almost_equal,
                        raises, assert_in)

from nose.plugins.skip import Skip, SkipTest
from .test_helpers import true_func, assert_equal_array_array, make_1d_traj

import numpy as np

import openpathsampling as paths
from openpathsampling.bias_function import *

import logging

quiet_loggers = ["initialization", "ensemble", "netcdfplus.objects",
                 "netcdfplus.netcdfplus", "pathmover", "netcdfplus.base"]
for logger in quiet_loggers:
    logging.getLogger("openpathsampling."+logger).setLevel(logging.CRITICAL)


class TestBiasEnsembleTable(object):
    def setup(self):
        # create the network
        paths.InterfaceSet._reset()
        xval = paths.FunctionCV(name="xA", f=lambda s : s.xyz[0][0])
        self.stateA = paths.CVDefinedVolume(xval, -1.0, -0.5).named("A")
        self.stateB = paths.CVDefinedVolume(xval, 0.5, float("inf")).named("B")
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
            self.change_12 : old_div(0.2, 0.5),
            self.change_10 : old_div(1.0, 0.5),
            self.change_20 : old_div(1.0, 0.2),
            self.change_21 : old_div(0.5, 0.2)
        }
        for change in list(change_vals.keys()):
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
            self.change_21 : old_div(0.2, 0.5),
            self.change_01 : old_div(1.0, 0.5),
            self.change_02 : old_div(1.0, 0.2),
            self.change_12 : old_div(0.5, 0.2)
        }
        for change in list(change_vals.keys()):
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
        ifacesA = self.ifacesA[:-1]
        xval2 = paths.FunctionCV(name="xB", f=lambda s : 0.5-s.xyz[0][0])
        ifacesB = paths.VolumeInterfaceSet(xval2, float("-inf"),
                                           [0.0, 0.1, 0.2])
        xval3 = paths.FunctionCV(name="xC", f=lambda s : s.xyz[0][0]-2.0)
        stateC = paths.CVDefinedVolume(self.xval, -3.0, 2.0)
        ifacesC = paths.VolumeInterfaceSet(xval3, -1.0, [0.0, 0.1, 0.2, 0.3])
        network = paths.MISTISNetwork(
            [(self.stateA, ifacesA, self.stateB),
             (self.stateB, ifacesB, self.stateA),
             (stateC, ifacesC, self.stateA)],
            ms_outers=paths.MSOuterTISInterface.from_lambdas(
                {ifacesA: -0.2, ifacesB: 0.3}
            )
        )
        ens_A = network.transitions[(self.stateA, self.stateB)].ensembles
        ens_B = network.transitions[(self.stateB, self.stateA)].ensembles
        ens_C = network.transitions[(stateC, self.stateA)].ensembles
        ms_outer = list(network.special_ensembles['ms_outer'].keys())[0]
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

class TestSRTISBiasFromNetwork(object):
    def setup(self):
        paths.InterfaceSet._reset()
        xval = paths.CoordinateFunctionCV(name="xA",
                                          f=lambda s : s.xyz[0][0])
        self.stateA = paths.CVDefinedVolume(xval, -1.0, -0.5).named("A")
        self.stateB = paths.CVDefinedVolume(xval, 0.5, float("inf")).named("B")
        self.ifacesA = paths.VolumeInterfaceSet(xval, -1.0, 
                                                [-0.5, -0.4, -0.3, -0.2])
        self.ifacesB = paths.VolumeInterfaceSet(xval, [0.5, 0.4, 0.3, 0.2],
                                                1.0)
        self.tcp_A = paths.numerics.LookupFunction(
            ordinate=[-0.5, -0.4, -0.3, -0.2, -0.1],
            abscissa=[1.0, 0.5, 0.25, 0.125, 0.0625]
        )
        self.tcp_B = paths.numerics.LookupFunction(
            ordinate=[0.5, 0.4, 0.3, 0.2, 0.1],
            abscissa=[1.0, 0.2, 0.04, 0.008, 0.0016]
        )


    def test_bias_from_network(self):
        # force the TCP in
        network = paths.MISTISNetwork([
            (self.stateA, self.ifacesA, self.stateB)
        ])
        network.sampling_transitions[0].tcp = self.tcp_A
        bias = paths.SRTISBiasFromNetwork(network)
        transition = list(network.transitions.values())[0]  # only one

        # check reciprocal of symmetric partners
        for i in range(4):
            for j in range(i, 4):
                assert_equal(bias.dataframe.loc[i, j],
                             old_div(1.0, bias.dataframe.loc[j, i]))

        for i in range(len(transition.ensembles) - 1):
            ens_to = transition.ensembles[i]
            ens_from = transition.ensembles[i + 1]
            assert_equal(bias.bias_value(ens_from, ens_to), 0.5)

        for i in range(len(transition.ensembles) - 2):
            ens_to = transition.ensembles[i]
            ens_from = transition.ensembles[i + 2]
            assert_equal(bias.bias_value(ens_from, ens_to), 0.25)

    @raises(RuntimeError)
    def test_fail_without_tcp(self):
        network = paths.MISTISNetwork([
            (self.stateA, self.ifacesA, self.stateB)
        ])
        bias = paths.SRTISBiasFromNetwork(network)

    @raises(RuntimeError)
    def test_fail_without_lambdas(self):
        fake_ifaceA = paths.InterfaceSet(cv=self.ifacesA.cv,
                                         volumes=self.ifacesA.volumes,
                                         direction=self.ifacesA.direction)
        network = paths.MISTISNetwork([
            (self.stateA, fake_ifaceA, self.stateB)
        ])
        network.sampling_transitions[0].tcp = self.tcp_A
        bias = paths.SRTISBiasFromNetwork(network)

    def test_bias_from_ms_network(self):
        ms_outer = paths.MSOuterTISInterface.from_lambdas(
            {self.ifacesA : -0.1, self.ifacesB : 0.1}
        )
        network = paths.MISTISNetwork(
            [(self.stateA, self.ifacesA, self.stateB),
             (self.stateB, self.ifacesB, self.stateA)],
            ms_outers=[ms_outer]
        )
        transition_AB = None
        transition_BA = None
        for t in network.sampling_transitions:
            if t.stateA == self.stateA:
                t.tcp = self.tcp_A
                transition_AB = t
            elif t.stateA == self.stateB:
                t.tcp = self.tcp_B
                transition_BA = t
            else:
                print([t.stateA, t.stateB])
                print([self.stateA, self.stateB])
                raise RuntimeError("Weird states in test transition")

        bias = paths.SRTISBiasFromNetwork(network)

        n_ensembles = len(bias.dataframe.index)
        for i in range(n_ensembles):
            for j in range(i, n_ensembles):
                if not np.isnan(bias.dataframe.loc[i, j]):
                    np.testing.assert_almost_equal(
                        bias.dataframe.loc[i, j],
                        old_div(1.0, bias.dataframe.loc[j, i])
                    )

        for i in range(len(transition_AB.ensembles) - 1):
            ens_to = transition_AB.ensembles[i]
            ens_from = transition_AB.ensembles[i + 1]
            assert_almost_equal(bias.bias_value(ens_from, ens_to), 0.5)

        for i in range(len(transition_BA.ensembles) - 1):
            ens_to = transition_BA.ensembles[i]
            ens_from = transition_BA.ensembles[i + 1]
            assert_almost_equal(bias.bias_value(ens_from, ens_to), 0.2)

        for ensA in transition_AB.ensembles:
            for ensB in transition_BA.ensembles:
                assert_equal(np.isnan(bias.bias_value(ensA, ensB)), True)
                assert_equal(np.isnan(bias.bias_value(ensB, ensA)), True)

        assert_almost_equal(bias.bias_value(transition_BA.ensembles[-1],
                                            network.ms_outers[0]),
                            old_div(5.0, 2))
        assert_almost_equal(bias.bias_value(transition_AB.ensembles[-1],
                                            network.ms_outers[0]),
                            old_div(2.0, 2))





        raise SkipTest
