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
                          make_1d_traj,
                          CalvinistDynamics
                         )

from openpathsampling.ensemble import LengthEnsemble
from openpathsampling.sample import SampleSet, Sample
from openpathsampling.pathmover import *

from openpathsampling.shooting import UniformSelector

from openpathsampling.volume import LambdaVolume
from test_helpers import CallIdentity
from openpathsampling.trajectory import Trajectory
from openpathsampling.snapshot import Snapshot
from openpathsampling.ensemble import EnsembleFactory as ef
from openpathsampling.orderparameter import OP_Function, OrderParameter

import logging
logging.getLogger('openpathsampling.pathmover').setLevel(logging.DEBUG)
logging.getLogger('openpathsampling.initialization').setLevel(logging.CRITICAL)

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

class testShootingMover(object):
    def setup(self):
        self.dyn = CalvinistDynamics([-0.1, 0.1, 0.3, 0.5, 0.7, 
                                      -0.1, 0.2, 0.4, 0.6, 0.8,
                                     ])
        PathMover.engine = self.dyn
        try:
            op = OP_Function("myid", fcn=lambda snap : 
                             Trajectory([snap])[0].coordinates()[0][0])
        except ValueError:
            op = OrderParameter.get_existing('myid')
        stateA = LambdaVolume(op, -100, 0.0)
        stateB = LambdaVolume(op, 0.65, 100)
        self.tps = ef.A2BEnsemble(stateA, stateB)
        init_traj = make_1d_traj(
            coordinates=[-0.1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7],
            velocities=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        )
        self.init_samp = SampleSet([Sample(
            trajectory=init_traj,
            replica=0,
            ensemble=self.tps
        )])

class testForwardShootMover(testShootingMover):
    def test_move(self):
        mover = ForwardShootMover(UniformSelector(), replicas=[0])
        self.dyn.initialized = True
        movepath = mover.move(self.init_samp)
        newsamp = self.init_samp + movepath
        assert_equal(len(newsamp), 1)
        assert_equal(newsamp[0].details.accepted, True)
        assert_equal(newsamp[0].ensemble(newsamp[0].trajectory), True)
        assert_equal(newsamp[0].trajectory, newsamp[0].details.trial)

class testBackwardShootMover(testShootingMover):
    def test_move(self):
        mover = BackwardShootMover(UniformSelector(), replicas=[0])
        self.dyn.initialized = True
        movepath = mover.move(self.init_samp)
        newsamp = self.init_samp + movepath
        assert_equal(len(newsamp), 1)
        assert_equal(newsamp[0].details.accepted, True)
        assert_equal(newsamp[0].ensemble(newsamp[0].trajectory), True)
        assert_equal(newsamp[0].trajectory, newsamp[0].details.trial)

class testOneWayShootingMover(testShootingMover):
    def test_mover_initialization(self):
        mover = OneWayShootingMover(UniformSelector, replicas=[0])
        assert_equal(len(mover.movers), 2)
        assert_equal(isinstance(mover, RandomChoiceMover), True)
        assert_equal(isinstance(mover, OneWayShootingMover), True)
        moverclasses = [m.__class__ for m in mover.movers]
        assert_equal(ForwardShootMover in moverclasses, True)
        assert_equal(BackwardShootMover in moverclasses, True)

class testPathReversalMover(object):
    def setup(self):
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
        samp = (gs_AXA + self.move.move(gs_AXA))[0]
        assert_equal(samp.details.accepted, True)

    def test_A_A_path(self):
        trajA_A = make_1d_traj(coordinates=[-0.3, 0.1, -0.4])
        sampA_A = Sample(trajectory=trajA_A,
                         ensemble=self.tis,
                         replica=0,
                         details=MoveDetails())
        gs_A_A = SampleSet([sampA_A])
        samp = (gs_A_A + self.move.move(gs_A_A))[0]
        assert_equal(samp.details.accepted, False)


    def test_AB_path(self):
        trajAXB = make_1d_traj(coordinates=[-0.2, 0.75, 1.8])
        sampAXB = Sample(trajectory=trajAXB,
                         ensemble=self.tis,
                         replica=0,
                         details=MoveDetails())
        gs_AXB = SampleSet([sampAXB])
        samp = (gs_AXB + self.move.move(gs_AXB))[0]
        assert_equal(samp.details.accepted, False)

    def test_BA_path(self):
        trajBXA = make_1d_traj(coordinates=[1.2, 0.7, -0.25])
        sampBXA = Sample(trajectory=trajBXA,
                         ensemble=self.tis,
                         replica=0,
                         details=MoveDetails())
        gs_BXA = SampleSet([sampBXA])
        samp = (gs_BXA + self.move.move(gs_BXA))[0]
        assert_equal(samp.details.accepted, True)

class testReplicaExchangeMover(object):
    def setup(self):
        try:
            op = OP_Function("myid", fcn=lambda snap : 
                             Trajectory([snap])[0].coordinates()[0][0])
        except ValueError:
            op = OrderParameter.get_existing('myid')
        state1 = LambdaVolume(op, -100, 0.0)
        state2 = LambdaVolume(op, 1, 100)
        volA = LambdaVolume(op, -100, 0.25)
        volB = LambdaVolume(op, -100, 0.50)
        self.tisA = ef.TISEnsemble(state1, state2, volA)
        self.tisB = ef.TISEnsemble(state1, state2, volB)
        self.traj0 = make_1d_traj([-0.1, 0.2, 0.3, 0.1, -0.2])
        self.traj1 = make_1d_traj([-0.1, 0.1, 0.4, 0.6, 0.3, 0.2, -0.15]) 
        self.traj2 = make_1d_traj([-0.1, 0.2, 0.3, 0.7, 0.6, 0.4, 0.1, -0.15])
        self.sampA0 = Sample(replica=0, trajectory=self.traj0, ensemble=self.tisA)
        self.sampB1 = Sample(replica=1, trajectory=self.traj1, ensemble=self.tisB)
        self.sampA2 = Sample(replica=2, trajectory=self.traj2, ensemble=self.tisA)
        self.gs_B1A2 = SampleSet([self.sampB1, self.sampA2])
        self.gs_A0B1 = SampleSet([self.sampA0, self.sampB1])

    def test_repex_ens_acc(self):
        repex_AB = ReplicaExchangeMover(ensembles=[[self.tisA, self.tisB]])
        samples_B2A1_ens = repex_AB.move(self.gs_B1A2)
        assert_equal(len(samples_B2A1_ens), 2)
        for sample in samples_B2A1_ens:
            assert_equal(sample.details.accepted, True)
            assert_equal(sample.trajectory, sample.details.result)
            assert_equal(sample.details.trial, sample.details.result)
        B2 = [s for s in samples_B2A1_ens if s.ensemble==self.tisB]
        assert_equal(len(B2), 1)
        assert_equal(B2[0].trajectory, self.traj2)
        A1 = [s for s in samples_B2A1_ens if s.ensemble==self.tisA]
        assert_equal(len(A1), 1)
        assert_equal(A1[0].trajectory, self.traj1)

    def test_repex_ens_rej(self):
        repex_AB = ReplicaExchangeMover(ensembles=[[self.tisA, self.tisB]])
        samples_A0B1_ens = repex_AB.move(self.gs_A0B1)
        assert_equal(len(samples_A0B1_ens), 2)
        for sample in samples_A0B1_ens:
            assert_equal(sample.details.accepted, False)
            assert_equal(sample.trajectory, sample.details.result)
            assert_not_equal(sample.details.trial, sample.details.result)
        A0 = [s for s in samples_A0B1_ens if s.ensemble==self.tisA]
        assert_equal(len(A0), 1)
        assert_equal(A0[0].trajectory, self.traj0)
        B1 = [s for s in samples_A0B1_ens if s.ensemble==self.tisB]
        assert_equal(len(B1), 1)
        assert_equal(B1[0].trajectory, self.traj1)


    def test_repex_rep_acc(self):
        repex_12 = ReplicaExchangeMover(replicas=[[1,2]])
        samples_B2A1_rep = repex_12.move(self.gs_B1A2)
        assert_equal(len(samples_B2A1_rep), 2)
        for sample in samples_B2A1_rep:
            assert_equal(sample.details.accepted, True)
            assert_equal(sample.trajectory, sample.details.result)
            assert_equal(sample.details.trial, sample.details.result)
        B2 = [s for s in samples_B2A1_rep if s.ensemble==self.tisB]
        assert_equal(len(B2), 1)
        assert_equal(B2[0].trajectory, self.traj2)
        A1 = [s for s in samples_B2A1_rep if s.ensemble==self.tisA]
        assert_equal(len(A1), 1)
        assert_equal(A1[0].trajectory, self.traj1)


class testRandomChoiceMover(object):
    def setup(self):
        traj = Trajectory([-0.5, 0.7, 1.1])
        op = CallIdentity()
        volA = LambdaVolume(op, -100, 0.0)
        volB = LambdaVolume(op, 1.0, 100)
        volX = LambdaVolume(op, -100, 0.25)
        self.tis = ef.TISEnsemble(volA, volB, volX)
        self.tps = ef.A2BEnsemble(volA, volB)
        self.len3 = LengthEnsemble(3)
        self.init_samp = SampleSet([Sample(trajectory=traj,
                                           ensemble=self.len3, 
                                           replica=0, 
                                           details=MoveDetails())])
        self.hop_to_tis = EnsembleHopMover(ensembles=[[self.len3, self.tis]])
        self.hop_to_tps = EnsembleHopMover(ensembles=[[self.len3, self.tps]])
        self.mover = RandomChoiceMover([self.hop_to_tis, self.hop_to_tps])

    def test_random_choice(self):
        # test that both get selected, but that we always return only one
        # sample
        count = {}
        for t in range(100):
            samples = self.mover.move(self.init_samp)
            assert_equal(len(samples), 1)
#            try:
                # Since self is the root mover, mover_path[-1] is self.
                # That means that mover_path[-2] is the mover that this
                # mover chose.
#                count[samples[0].details.mover_path[-2]] += 1
#            except KeyError:
#                count[samples[0].details.mover_path[-2]] = 1
#        assert_equal(len(count.keys()), 2)

    def test_restricted_by_replica(self):
        raise SkipTest

    def test_restricted_by_ensemble(self):
        raise SkipTest

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
        movepath = move.move(gs)
        samples = movepath.changes
        assert_equal(len(samples), 3)
        for sample in samples:
            assert_equal(sample.details.accepted, True)
        gs = gs + movepath
        assert_equal(gs[0].ensemble, self.tps)

    def test_first_rejected(self):
        move = SequentialMover(movers=self.first_rejected_movers)
        gs = SampleSet(self.init_sample)
        movepath = move.move(gs)
        samples = movepath.changes
        # @DWHS: This should have two samples since two are accepted
        # and thus applied
        assert_equal(len(samples), 2)

        allsamp = movepath.all_samples
        assert_equal(allsamp[0].details.accepted, False)
        assert_equal(allsamp[1].details.accepted, True)
        assert_equal(allsamp[2].details.accepted, True)
        gs = gs + movepath
        assert_equal(gs[0].ensemble, self.tps)

    def test_last_rejected(self):
        move = SequentialMover(movers=self.last_rejected_movers)
        gs = SampleSet(self.init_sample)
        movepath = move.move(gs)
        samples = movepath.changes
        assert_equal(len(samples), 2)
        # @DWHS: I think if the last is rejected then there should only be two
        # samples to be used, since the last one is not accepted and thus
        # discarded (does not mean that it is not stored!!!)
        allsamp = movepath.all_samples
        assert_equal(allsamp[0].details.accepted, True)
        assert_equal(allsamp[1].details.accepted, True)
        assert_equal(allsamp[2].details.accepted, False)
        gs = gs + movepath
        assert_equal(gs[0].ensemble, self.tps)

    def test_restricted_by_replica(self):
        raise SkipTest

    def test_restricted_by_ensemble(self):
        raise SkipTest

class testPartialAcceptanceSequentialMover(testSequentialMover):
    def test_everything_accepted(self):
        move = PartialAcceptanceSequentialMover(movers=self.everything_accepted_movers)
        gs = SampleSet(self.init_sample)
        movepath = move.move(gs)
        samples = movepath.changes
        assert_equal(len(samples), 3)
        for sample in samples:
            assert_equal(sample.details.accepted, True)
        assert_equal(len(movepath.all_samples,),3)
        gs = gs + movepath
        assert_equal(gs[0].ensemble, self.tps)

    def test_first_rejected(self):
        move = PartialAcceptanceSequentialMover(movers=self.first_rejected_movers)
        gs = SampleSet(self.init_sample)
        movepath = move.move(gs)
        samples = movepath.changes
        # returns zero sample since even the first is rejected
        # the first one is still stored
        assert_equal(len(samples), 0)
        allsamp = movepath.all_samples
        assert_equal(len(allsamp), 1)
        assert_equal(allsamp[0].details.accepted, False)
        gs = gs + movepath
        assert_equal(gs[0].ensemble, self.len3)

    def test_last_rejected(self):
        move = PartialAcceptanceSequentialMover(movers=self.last_rejected_movers)
        gs = SampleSet(self.init_sample)
        movepath = move.move(gs)
        samples = movepath.changes
        # @see above, this should return 2 samples. Important the third is
        # still run!
        assert_equal(len(samples), 2)
        allsamp = movepath.all_samples
        assert_equal(len(allsamp), 3)

        assert_equal(allsamp[0].details.accepted, True)
        assert_equal(allsamp[1].details.accepted, True)
        assert_equal(allsamp[2].details.accepted, False)
        gs = gs + movepath
        assert_equal(gs[0].ensemble, self.tps)

    def test_restricted_by_replica(self):
        raise SkipTest

    def test_restricted_by_ensemble(self):
        raise SkipTest

class testConditionalSequentialMover(testSequentialMover):
    def test_everything_accepted(self):
        move = ConditionalSequentialMover(movers=self.everything_accepted_movers)
        gs = SampleSet(self.init_sample)
        movepath = move.move(gs)
        samples = movepath.changes
        assert_equal(len(samples), 3)
        for sample in samples:
            assert_equal(sample.details.accepted, True)
        gs = gs + movepath
        assert_equal(gs[0].ensemble, self.tps)

    def test_first_rejected(self):
        move = ConditionalSequentialMover(movers=self.first_rejected_movers)
        gs = SampleSet(self.init_sample)
        movepath = move.move(gs)
        samples = movepath.changes
        # should be zero since the move is completely rejected
        assert_equal(len(samples), 0)
        allsamp = movepath.all_samples
        assert_equal(len(allsamp), 1)
        assert_equal(allsamp[0].details.accepted, False)
        gs = gs + movepath
        assert_equal(gs[0].ensemble, self.len3)

    def test_last_rejected(self):
        move = ConditionalSequentialMover(movers=self.last_rejected_movers)
        gs = SampleSet(self.init_sample)
        movepath = move.move(gs)
        samples = movepath.changes
        # number of accepted samples is 0 for this type of mover
        assert_equal(len(samples), 0)
        allsamp = movepath.all_samples
        assert_equal(len(allsamp), 3)

        # check here if last actual samples was false
        # this actually allows to see later if the single samples were
        # accepted or not, even from the movepath without loading samples
        assert_equal(allsamp[0].details.accepted, True)
        assert_equal(allsamp[1].details.accepted, True)
        assert_equal(allsamp[2].details.accepted, False)
        gs = gs + movepath
        assert_equal(gs[0].ensemble, self.len3)

    def test_restricted_by_replica(self):
        raise SkipTest

    def test_restricted_by_ensemble(self):
        raise SkipTest

