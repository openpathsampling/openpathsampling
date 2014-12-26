'''
@author: David W.H. Swenson
'''

import os
import numpy as np

from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        assert_almost_equal, raises)
from nose.plugins.skip import Skip, SkipTest
from test_helpers import (assert_equal_array_array, 
                          assert_not_equal_array_array)

from opentis.ensemble import LengthEnsemble
from opentis.sample import SampleSet, Sample
from opentis.pathmover import *

from opentis.volume import LambdaVolume
from test_helpers import CallIdentity
from opentis.ensemble import EnsembleFactory as ef

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
        pass

    def test_AA_path(self):
        pass

    def test_AB_path(self):
        pass

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
        traj = [-0.5, 0.7, 1.1]
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


    def test_everything_accepted(self):
        move = SequentialMover(movers=[self.hop_to_tis,
                                       self.hop_to_len3,
                                       self.hop_to_tps])
        gs = GlobalState(self.init_sample)
        samples = move.move(gs)
        assert_equal(len(samples), 3)
        for sample in samples:
            assert_equal(sample.details.accepted, True)

        

    def test_first_rejected(self):
        pass

    def test_last_rejected(self):
        pass

    def test_restricted_by_replica(self):
        pass

    def test_restricted_by_ensemble(self):
        pass

class testPartialAcceptanceSequentialMover(testSequentialMover):
    def test_everything_accepted(self):
        pass

    def test_first_rejected(self):
        pass

    def test_last_rejected(self):
        pass

    def test_restricted_by_replica(self):
        pass

    def test_restricted_by_ensemble(self):
        pass

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
        pass

    def test_restricted_by_ensemble(self):
        pass

