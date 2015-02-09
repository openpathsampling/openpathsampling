'''
@author: David W.H. Swenson
'''

from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        assert_almost_equal, raises)
from nose.plugins.skip import Skip, SkipTest
from test_helpers import (assert_equal_array_array, items_equal,
                          assert_not_equal_array_array,
                          make_1d_traj,
                          CalvinistDynamics,
                          CallIdentity
                         )
import openpathsampling as paths
from openpathsampling.ensemble import LengthEnsemble
from openpathsampling.pathmover import *

from openpathsampling.shooting import UniformSelector

from openpathsampling.volume import LambdaVolume
from test_helpers import CallIdentity
from openpathsampling.trajectory import Trajectory
from openpathsampling.ensemble import EnsembleFactory as ef
from openpathsampling.orderparameter import OP_Function, OrderParameter

import logging
#logging.getLogger('openpathsampling.pathmover').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.ensemble').setLevel(logging.CRITICAL)


#logging.getLogger('openpathsampling.pathmover').propagate = False
#logging.getLogger('openpathsampling.initialization').propagate = False

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

def assert_sampleset_accepted(sampleset, results):
    for sample, result in zip(sampleset, results):
        assert_equal(sample.details.accepted, result)


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
                             snap.coordinates[0][0])
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
                             snap.coordinates[0,0])
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

class testReplicaIDChangeMover(object):
    def setup(self):
        pass

    def test_replica_in_sampleset(self):
        raise SkipTest

    def test_replica_not_in_sampleset(self):
        raise SkipTest


class testReplicaExchangeMover(object):
    def setup(self):
        try:
            op = OP_Function("myid", fcn=lambda snap : 
                             snap.coordinates[0][0])
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
        repex_movepath = repex_AB.move(self.gs_A0B1)

        assert_equal(len(repex_movepath.samples), 0) # since rejected

        samples_A0B1_ens = repex_movepath.all_samples
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
        samples = movepath.samples
        assert_equal(len(samples), 3)
        for sample in samples:
            assert_equal(sample.details.accepted, True)
        gs = gs + movepath
        assert_equal(gs[0].ensemble, self.tps)

    def test_first_rejected(self):
        move = SequentialMover(movers=self.first_rejected_movers)
        gs = SampleSet(self.init_sample)
        movepath = move.move(gs)
        samples = movepath.samples
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
        samples = movepath.samples
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
        samples = movepath.samples
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
        samples = movepath.samples
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
        samples = movepath.samples
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
        samples = movepath.samples
        assert_equal(len(samples), 3)
        for sample in samples:
            assert_equal(sample.details.accepted, True)
        gs = gs + movepath
        assert_equal(gs[0].ensemble, self.tps)

    def test_first_rejected(self):
        move = ConditionalSequentialMover(movers=self.first_rejected_movers)
        gs = SampleSet(self.init_sample)
        movepath = move.move(gs)
        samples = movepath.samples
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
        samples = movepath.samples
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

class SubtrajectorySelectTester(object):

    def setup(self):
        op = CallIdentity()
        vol = paths.LambdaVolume(op, -0.5, 0.5)
        inX = paths.InXEnsemble(vol)
        outX = paths.OutXEnsemble(vol)
        self.ensemble = paths.SequentialEnsemble([
            inX, outX, inX, outX, inX, outX, inX
        ])
        self.subensemble = paths.SequentialEnsemble([
            paths.SingleFrameEnsemble(inX),
            outX,
            paths.SingleFrameEnsemble(inX)
        ])
        self.traj_with_3_subtrajs = Trajectory(
            [0.0, 0.0, 1.0, 1.0, 0.0, 1.0, 0.0, 0.0, 2.0, 0.0]
        )
        self.subtraj0 = Trajectory([0.0, 1.0, 1.0, 0.0])
        self.subtraj1 = Trajectory([0.0, 1.0, 0.0])
        self.subtraj2 = Trajectory([0.0, 2.0, 0.0])
        self.gs = SampleSet(Sample(
            replica=0,
            ensemble=self.subensemble,
            trajectory=self.traj_with_3_subtrajs
        ))

    def test_paths_in_ensemble(self):
        # more a test of SequentialEnsemble, but also a test of sanity
        # before the real tests
        assert_equal(self.ensemble(self.traj_with_3_subtrajs), True)
        assert_equal(self.subensemble(self.subtraj0), True)
        assert_equal(self.subensemble(self.subtraj1), True)
        assert_equal(self.subensemble(self.subtraj2), True)

class testRandomSubtrajectorySelectMover(SubtrajectorySelectTester):
    def test_accepts_all(self):
        mover = RandomSubtrajectorySelectMover(self.subensemble)
        found = {}
        for t in range(100):
            movepath = mover.move(self.gs)
            samples = movepath.samples
            assert_equal(len(samples), 1)
            assert_equal(self.subensemble, samples[0].ensemble)
            assert_equal(self.subensemble(samples[0].trajectory), True)
            assert_equal(self.ensemble(samples[0].trajectory), False)
            if samples[0].trajectory == self.subtraj0:
                found[0] = True
            elif samples[0].trajectory == self.subtraj1:
                found[1] = True
            elif samples[0].trajectory == self.subtraj2:
                found[2] = True
            else:
                raise RuntimeError("Subtraj unknown!")
        assert_equal(found[0] and found[1] and found[2], True)

    def test_nl_fails(self):
        raise SkipTest

    def test_nothing_allowed(self):
        mover = RandomSubtrajectorySelectMover(self.subensemble)
        traj_with_no_subtrajs = Trajectory([0.0, 0.0, 0.0])
        self.gs[0].trajectory = traj_with_no_subtrajs
        movepath = mover.move(self.gs)
        samples = movepath.samples
        assert_equal(samples[0].trajectory, paths.Trajectory([]))

class testFirstSubtrajectorySelectMover(SubtrajectorySelectTester):
    def test_move(self):
        mover = FirstSubtrajectorySelectMover(self.subensemble)
        movepath = mover.move(self.gs)
        samples = movepath.samples
        assert_equal(len(samples), 1)
        assert_equal(self.subensemble, samples[0].ensemble)
        assert_equal(self.subensemble(samples[0].trajectory), True)
        assert_equal(self.ensemble(samples[0].trajectory), False)
        assert_equal(samples[0].trajectory, self.subtraj0)

class testFinalSubtrajectorySelectMover(SubtrajectorySelectTester):
    def test_move(self):
        mover = FinalSubtrajectorySelectMover(self.subensemble)
        movepath = mover.move(self.gs)
        samples = movepath.samples
        assert_equal(len(samples), 1)
        assert_equal(self.subensemble, samples[0].ensemble)
        assert_equal(self.subensemble(samples[0].trajectory), True)
        assert_equal(self.ensemble(samples[0].trajectory), False)
        assert_equal(samples[0].trajectory, self.subtraj2)

class testForceEnsembleChangeMover(object):
    def setup(self):
        traj = Trajectory([-0.5, 0.7, 1.1])
        op = CallIdentity()
        volA = LambdaVolume(op, -100, 0.0)
        volB = LambdaVolume(op, 1.0, 100)
        volX = LambdaVolume(op, -100, 0.25)
        self.tis = ef.TISEnsemble(volA, volB, volX)
        self.len3 = LengthEnsemble(3)
        self.len2 = LengthEnsemble(2)
        self.gs = SampleSet(Sample(
            trajectory=traj,
            ensemble=self.tis,
            replica=0
        ))

    def test_in_ensemble(self):
        mover = ForceEnsembleChangeMover(ensembles=[[self.tis, self.len3]])
        movepath = mover.move(self.gs)
        samples = movepath.samples
        assert_equal(samples[0].details.initial_ensemble(samples[0].trajectory),
                     True)
        assert_equal(samples[0].ensemble(samples[0].trajectory), True)
        assert_equal(samples[0].ensemble, self.len3)

    def test_not_in_ensemble(self):
        mover = ForceEnsembleChangeMover(ensembles=[[self.tis, self.len2]])
        movepath = mover.move(self.gs)
        samples = movepath.samples
        assert_equal(samples[0].details.initial_ensemble(samples[0].trajectory),
                     True)
        assert_equal(samples[0].ensemble, self.len2)
        assert_equal(samples[0].ensemble(samples[0].trajectory), False)

class testMinusMover(object):
    def setup(self):
        try:
            op = OP_Function("myid", fcn=lambda snap : snap.coordinates[0,0])
        except ValueError:
            op = OrderParameter.get_existing('myid')
        volA = LambdaVolume(op, -100, 0.0)
        volB = LambdaVolume(op, 1.0, 100)
        volX = LambdaVolume(op, -100, 0.25)
        self.dyn = CalvinistDynamics([
            # successful move: (backward extension then forward)
            -0.13, 0.13, 0.33, -0.11, -0.12, 0.12, 0.32, -0.131,
            # never leaves state: 
            -0.15, -0.15, -0.15, -0.15, -0.15, -0.15, -0.15, -0.15, -0.15,
            -0.15, -0.15, -0.15, -0.15, -0.15, -0.15, -0.15, -0.15, -0.15,
            -0.25, 
            -0.15, -0.15, -0.15, -0.15, -0.15, -0.15, -0.15, -0.15, -0.15,
            -0.15, -0.15, -0.15, -0.15, -0.15, -0.15, -0.15, -0.15, -0.15,
            # goes to other state:
            1.16, 1.26, 1.16, -0.16, 1.16, 1.26, 1.16
        ])
        PathMover.engine = self.dyn
        self.dyn.initialized = True
        self.innermost = ef.TISEnsemble(volA, volB, volX)
        self.minus = paths.MinusInterfaceEnsemble(volA, volX)
        self.mover = MinusMover(minus_ensemble=self.minus,
                                innermost_ensemble=self.innermost)
        self.first_segment = [-0.1, 0.1, 0.3, 0.1, -0.15] 
        self.list_innermost = [-0.11, 0.11, 0.31, 0.11, -0.12]
        self.second_segment = [-0.25, 0.2, 0.4, 0.2, -0.2]
        init_minus = make_1d_traj(
            coordinates=self.first_segment + [-0.35] + self.second_segment,
            velocities=[1.0]*11
        )
        self.minus_sample = Sample(
            replica=-1,
            trajectory=init_minus,
            ensemble=self.minus
        )

    def test_setup_sanity(self):
        # sanity checks to make sure that what we set up makes sense
        assert_equal(self.minus_sample.ensemble(self.minus_sample.trajectory),
                    True)
        first_subtraj = FirstSubtrajectorySelectMover(
            subensemble=self.minus._segment_ensemble
        )
        movepath = first_subtraj.move(SampleSet(self.minus_sample))
        samples = movepath.samples
        assert_equal(samples[0].ensemble(samples[0].trajectory), True)
        final_subtraj = FinalSubtrajectorySelectMover(
            subensemble=self.minus._segment_ensemble
        )
        movepath = final_subtraj.move(SampleSet(self.minus_sample))
        samples = movepath.samples
        assert_equal(samples[0].ensemble(samples[0].trajectory), True)
        assert_equal(samples[0].ensemble, self.minus._segment_ensemble)
        

    def test_successful_move(self):
        init_innermost = make_1d_traj(self.list_innermost, [1.0]*5)
        init_sample = Sample(
            replica=0,
            trajectory=init_innermost,
            ensemble=self.innermost
        )
        gs = SampleSet([init_sample, self.minus_sample])

        extend_forward =  self.list_innermost + [0.12, 0.32, -0.131]
        extend_backward = [-0.13, 0.13, 0.33] + self.list_innermost

        assert_equal(self.minus(make_1d_traj(extend_forward)), True)
        assert_equal(self.minus(make_1d_traj(extend_backward)), True)

        seg_dir = {}
        for i in range(100):
            movepath = self.mover.move(gs).opened
            samples = movepath.samples
            assert_equal(len(samples), 5)
            s_inner = [s for s in samples if s.ensemble==self.innermost]
            s_minus = [s for s in samples if s.ensemble==self.minus]
            s_sub = [s for s in samples if s.ensemble==self.minus._segment_ensemble]
            assert_equal(len(s_inner), 1)
            assert_equal(len(s_minus), 2)
            assert_equal(len(s_sub), 2)

            for s in samples:
                assert_equal(s.details.accepted, True)

            key = ""
            s_inner0_xvals = [s.coordinates[0,0] for s in s_inner[0].trajectory]
            if items_equal(s_inner0_xvals, self.first_segment):
                key += "1"
            elif items_equal(s_inner0_xvals, self.second_segment):
                key += "2"
            else:
                print s_inner0_xvals
                raise RuntimeError("Chosen segment neither first nor last!")

            # final sample s_minus is accepted
            s_minus_xvals = [s.coordinates[0,0] for s in s_minus[-1].trajectory]
            if items_equal(s_minus_xvals, extend_forward):
                key += "f"
            elif items_equal(s_minus_xvals, extend_backward):
                key += "b"
            else:
                print s_minus_xvals
                raise RuntimeError("Unexpected minus extension result!")

            try:
                seg_dir[key] += 1
            except KeyError:
                seg_dir[key] = 1
        assert_equal(len(seg_dir.keys()), 4)

    def test_repex_fails_other_ensemble(self):
        innermost_other_ensemble = make_1d_traj([-0.11, 0.1, -0.12])
        samp_other_ensemble = Sample(
            replica=0,
            trajectory=innermost_other_ensemble,
            ensemble=self.innermost
        )
        gs = SampleSet([samp_other_ensemble, self.minus_sample])
        
        movepath = self.mover.move(gs).opened
        samples = movepath.all_samples
        assert_equal(self.innermost(innermost_other_ensemble), False)
        assert_equal(len(samples), 3) # stop after failed repex
        assert_equal(samples[0].details.accepted, True)
        assert_equal(samples[1].details.accepted, False)
        assert_equal(samples[2].details.accepted, False)

    def test_repex_fails_innermost_crosses_state(self):
        innermost_crosses_to_state = make_1d_traj([-0.11, 0.5, 1.8])
        samp_crosses_to_state = Sample(
            replica=0,
            trajectory=innermost_crosses_to_state,
            ensemble=self.innermost
        )
        gs = SampleSet([samp_crosses_to_state, self.minus_sample])
        
        movepath = self.mover.move(gs).opened
        samples = movepath.all_samples
        assert_equal(self.innermost(innermost_crosses_to_state), True)
        assert_equal(len(samples), 3) # stop after failed repex
        assert_sampleset_accepted(samples, [True, False, False])

    def test_repex_fails_minus_crosses_to_state(self):
        minus_crosses_to_state = make_1d_traj(
            [-0.11, 0.5, 1.8, 0.6, -0.12, 0.7, 1.7, 0.4, -0.13]
        )
        badminus_sample = Sample(
            replica=-1,
            trajectory=minus_crosses_to_state,
            ensemble=self.minus
        )
        init_sample = Sample(
            replica=0,
            trajectory=make_1d_traj(self.list_innermost, [1.0]*5),
            ensemble=self.innermost
        )
        gs = SampleSet([badminus_sample, init_sample])

        assert_equal(self.minus(minus_crosses_to_state), True)

        movepath = self.mover.move(gs).opened
        samples = movepath.all_samples
        assert_equal(len(samples), 3) # stop after failed repex
        assert_sampleset_accepted(samples, [True, False, False])

    def test_extension_fails(self):
        innermost_bad_extension = [-0.25, 0.1, 0.5, 0.1, -0.25]
        traj_bad_extension = make_1d_traj(innermost_bad_extension, [1.0]*5)
        samp_bad_extension = Sample(
            replica=0,
            trajectory=traj_bad_extension,
            ensemble=self.innermost
        )
        
        assert_equal(self.innermost(traj_bad_extension), True)

        gs = SampleSet([self.minus_sample, samp_bad_extension])
        movepath = self.mover.move(gs).opened
        samples = movepath.all_samples
        assert_equal(len(samples), 5) # reject the last one
        assert_sampleset_accepted(samples, [True] * 4 + [False])
        # this only happens due to length
        assert_equal(len(samples[-1].details.trial),
                     len(traj_bad_extension)+self.dyn.n_frames_max-1)
