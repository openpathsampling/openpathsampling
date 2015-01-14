'''
@author: David W.H. Swenson
'''

import os
import numpy as np

from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        assert_almost_equal, raises)
from nose.plugins.skip import Skip, SkipTest
from test_helpers import (assert_equal_array_array, 
                          assert_not_equal_array_array,
                          make_1d_traj
                         )

from opentis.ensemble import LengthEnsemble
from opentis.sample import SampleSet, Sample
from opentis.pathmover import *

from opentis.volume import LambdaVolume
from test_helpers import CallIdentity
from opentis.trajectory import Trajectory
from opentis.snapshot import Snapshot
from opentis.ensemble import EnsembleFactory as ef
from opentis.orderparameter import OP_Function, OrderParameter

import logging
logging.getLogger('opentis.pathmover').setLevel(logging.DEBUG)
logging.getLogger('opentis.initialization').setLevel(logging.CRITICAL)

class testMakeListOfPairs(object):
    def setup(self):
        self.correct = [ [0, 1], [2, 3], [4, 5] ]

    @raises(TypeError)
    def test_not_iterable_type_error(self):
        result = make_list_of_pairs(1)
    
    def test_list_of_list_pairs(self):
        result = make_list_of_pairs([[0, 1], [2, 3], [4, 5]])
        assert_equal_array_array(result, self.correct)

    @raises(AssertionError)
    def test_list_of_list_notpairs(self):
        result = make_list_of_pairs([[0], [1], [2]])

    @raises(AssertionError)
    def test_list_not_even(self):
        result = make_list_of_pairs([0, 1, 2, 3, 4])

    def test_list_even(self):
        result = make_list_of_pairs([0, 1, 2, 3, 4, 5])
        assert_equal_array_array(result, self.correct)

    def test_empty(self):
        assert_equal(make_list_of_pairs(None), None)

class testPathMover(object):
    def setup(self):
        self.l1 = LengthEnsemble(1)
        self.l2 = LengthEnsemble(2)
        self.l3 = LengthEnsemble(3)
        self.repsAll_ensNone = PathMover(replicas='all')
        self.reps12_ensNone = PathMover(replicas=[1, 2])
        self.repsAll_ens1 = PathMover(ensembles=self.l1)
        self.repsAll_ens12 = PathMover(ensembles=[self.l1, self.l2])
        self.reps1_ens2 = PathMover(replicas=1, ensembles=[self.l2])
        self.s1 = Sample(replica=1, ensemble=self.l2)
        self.s2 = Sample(replica=2, ensemble=self.l1)
        self.s3 = Sample(replica=3, ensemble=self.l1)
        self.s4 = Sample(replica=2, ensemble=self.l3)
        self.sset = SampleSet([self.s1, self.s2, self.s3, self.s4])

    def test_legal_sample_set(self):
        assert_items_equal(self.repsAll_ensNone.legal_sample_set(self.sset),
                           [self.s1, self.s2, self.s3, self.s4])
        assert_items_equal(self.reps12_ensNone.legal_sample_set(self.sset),
                           [self.s1, self.s2, self.s4])
        assert_items_equal(self.repsAll_ens12.legal_sample_set(self.sset),
                           [self.s1, self.s2, self.s3])
        assert_items_equal(self.repsAll_ens1.legal_sample_set(self.sset),
                           [self.s2, self.s3])
        assert_items_equal(self.reps1_ens2.legal_sample_set(self.sset),
                           [self.s1])
        assert_items_equal(
            self.repsAll_ensNone.legal_sample_set(self.sset, ensembles=self.l1),
            [self.s2, self.s3]
        )
        assert_items_equal(
            self.repsAll_ensNone.legal_sample_set(self.sset, ensembles=[self.l1]),
            [self.s2, self.s3]
        )


    def test_select_sample(self):
        assert_equal(self.reps1_ens2.select_sample(self.sset), self.s1)
        selected = self.repsAll_ens1.select_sample(self.sset)
        try:
            assert_equal(selected, self.s2)
        except AssertionError:
            assert_equal(selected, self.s3)

class testForwardShootMover(object):
    def setup(self):
        pass

    def test_move(self):
        raise SkipTest

class testBackwardShootMover(object):
    def setup(self):
        pass

    def test_move(self):
        raise SkipTest

class testPathReversalMover(object):
    def setup(self):
        #op = OrderParameter()
        try:
            op = OP_Function("myid", fcn=lambda snap : 
                             Trajectory([snap])[0].coordinates()[0][0])
        except ValueError:
            op = OrderParameter.get_existing('myid')
        volA = LambdaVolume(op, -100, 0.0)
        volB = LambdaVolume(op, 1.0, 100)
        volX = LambdaVolume(op, -100, 0.25)
        self.tis = ef.TISEnsemble(volA, volB, volX)
        self.move = PathReversalMover()
        self.op = op

    def test_AXA_path(self):
        trajAXA = make_1d_traj(coordinates=[-0.1, 0.75, -0.6],
                               velocities=[0.1, 0.05, -0.05])
        sampAXA = Sample(trajectory=trajAXA,
                         ensemble=self.tis,
                         replica=0,
                         details=MoveDetails())
        gs_AXA = SampleSet([sampAXA])
        samp = self.move.move(gs_AXA)[0]
        assert_equal(samp.details.accepted, True)

    def test_A_A_path(self):
        trajA_A = make_1d_traj(coordinates=[-0.3, 0.1, -0.4])
        sampA_A = Sample(trajectory=trajA_A,
                         ensemble=self.tis,
                         replica=0,
                         details=MoveDetails())
        gs_A_A = SampleSet([sampA_A])
        samp = self.move.move(gs_A_A)[0]
        assert_equal(samp.details.accepted, False)


    def test_AB_path(self):
        trajAXB = make_1d_traj(coordinates=[-0.2, 0.75, 1.8])
        sampAXB = Sample(trajectory=trajAXB,
                         ensemble=self.tis,
                         replica=0,
                         details=MoveDetails())
        gs_AXB = SampleSet([sampAXB])
        samp = self.move.move(gs_AXB)[0]
        assert_equal(samp.details.accepted, False)

    def test_BA_path(self):
        trajBXA = make_1d_traj(coordinates=[1.2, 0.7, -0.25])
        sampBXA = Sample(trajectory=trajBXA,
                         ensemble=self.tis,
                         replica=0,
                         details=MoveDetails())
        gs_BXA = SampleSet([sampBXA])
        samp = self.move.move(gs_BXA)[0]
        assert_equal(samp.details.accepted, True)


class testMixedMover(object):
    def setup(self):
        pass

    def test_both_get_selected(self):
        pass

    def test_only_one_gets_run(self):
        pass

    def test_restricted_by_replica(self):
        pass

    def test_restricted_by_ensemble(self):
        pass

class testSequentialMover(object):
    def setup(self):
        traj = Trajectory([-0.5, 0.7, 1.1])
        op = CallIdentity()
        volA = LambdaVolume(op, -100, 0.0)
        volB = LambdaVolume(op, 1.0, 100)
        volX = LambdaVolume(op, -100, 0.25)
        tis = ef.TISEnsemble(volA, volB, volX)
        tps = ef.A2BEnsemble(volA, volB)
        len3 = LengthEnsemble(3)
        len2 = LengthEnsemble(2)
        self.hop_to_tis = EnsembleHopMover(ensembles=[[tis, tis],
                                                      [tps, tis],
                                                      [len3, tis],
                                                      [len2, tis]])
        self.hop_to_tps = EnsembleHopMover(ensembles=[[tis, tps],
                                                      [tps, tps],
                                                      [len3, tps],
                                                      [len2, tps]])
        self.hop_to_len3 = EnsembleHopMover(ensembles=[[tis, len3],
                                                       [tps, len3],
                                                       [len3, len3],
                                                       [len2, len3]]) 
        self.hop_to_len2 = EnsembleHopMover(ensembles=[[tis, len2],
                                                      [tps, len2],
                                                      [len3, len2],
                                                      [len2, len2]])
        self.init_sample = Sample(trajectory=traj,
                                  ensemble=len3,
                                  replica=0,
                                  details=MoveDetails())
        self.tis = tis
        self.tps = tps
        self.len3 = len3
        self.len2 = len2
        self.everything_accepted_movers = [
            self.hop_to_tis, self.hop_to_len3, self.hop_to_tps
        ]
        self.first_rejected_movers = [
            self.hop_to_len2, self.hop_to_len3, self.hop_to_tps
        ]
        self.last_rejected_movers = [
            self.hop_to_tis, self.hop_to_tps, self.hop_to_len2
        ]

    def test_everything_accepted(self):
        move = SequentialMover(movers=self.everything_accepted_movers)
        gs = SampleSet(self.init_sample)
        samples = move.move(gs)
        assert_equal(len(samples), 3)
        for sample in samples:
            assert_equal(sample.details.accepted, True)
        gs.apply_samples(samples)
        assert_equal(gs[0].ensemble, self.tps)

    def test_first_rejected(self):
        move = SequentialMover(movers=self.first_rejected_movers)
        gs = SampleSet(self.init_sample)
        samples = move.move(gs)
        assert_equal(len(samples), 3)
        assert_equal(samples[0].details.accepted, False)
        assert_equal(samples[1].details.accepted, True)
        assert_equal(samples[2].details.accepted, True)
        gs.apply_samples(samples)
        assert_equal(gs[0].ensemble, self.tps)

    def test_last_rejected(self):
        move = SequentialMover(movers=self.last_rejected_movers)
        gs = SampleSet(self.init_sample)
        samples = move.move(gs)
        assert_equal(len(samples), 3)
        assert_equal(samples[0].details.accepted, True)
        assert_equal(samples[1].details.accepted, True)
        assert_equal(samples[2].details.accepted, False)
        gs.apply_samples(samples)
        assert_equal(gs[0].ensemble, self.tps)

    def test_restricted_by_replica(self):
        raise SkipTest

    def test_restricted_by_ensemble(self):
        raise SkipTest

class testPartialAcceptanceSequentialMover(testSequentialMover):
    def test_everything_accepted(self):
        move = PartialAcceptanceSequentialMover(movers=self.everything_accepted_movers)
        gs = SampleSet(self.init_sample)
        samples = move.move(gs)
        assert_equal(len(samples), 3)
        for sample in samples:
            assert_equal(sample.details.accepted, True)
        gs.apply_samples(samples)
        assert_equal(gs[0].ensemble, self.tps)

    def test_first_rejected(self):
        move = PartialAcceptanceSequentialMover(movers=self.first_rejected_movers)
        gs = SampleSet(self.init_sample)
        samples = move.move(gs)
        assert_equal(len(samples), 1)
        assert_equal(samples[0].details.accepted, False)
        gs.apply_samples(samples)
        assert_equal(gs[0].ensemble, self.len3)

    def test_last_rejected(self):
        move = PartialAcceptanceSequentialMover(movers=self.last_rejected_movers)
        gs = SampleSet(self.init_sample)
        samples = move.move(gs)
        assert_equal(len(samples), 3)
        assert_equal(samples[0].details.accepted, True)
        assert_equal(samples[1].details.accepted, True)
        assert_equal(samples[2].details.accepted, False)
        gs.apply_samples(samples)
        assert_equal(gs[0].ensemble, self.tps)

    def test_restricted_by_replica(self):
        raise SkipTest

    def test_restricted_by_ensemble(self):
        raise SkipTest

class testConditionalSequentialMover(testSequentialMover):
    def setup(self):
        pass

    def test_everything_accepted(self):
        pass

    def test_first_rejected(self):
        pass

    def test_last_rejected(self):
        pass

    def test_restricted_by_replica(self):
        raise SkipTest

    def test_restricted_by_ensemble(self):
        raise SkipTest

