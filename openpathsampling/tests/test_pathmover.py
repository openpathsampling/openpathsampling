'''
@author: David W.H. Swenson
'''
from __future__ import print_function
from __future__ import absolute_import

from builtins import zip
from builtins import str
from builtins import range
from builtins import object
from nose.plugins.skip import SkipTest
from nose.tools import (assert_equal, assert_not_equal, raises, assert_true,
                        assert_in, assert_not_in)
from numpy.testing import assert_allclose

from openpathsampling.collectivevariable import FunctionCV
from openpathsampling.engines.trajectory import Trajectory
from openpathsampling.ensemble import LengthEnsemble
from openpathsampling.pathmover import *
from openpathsampling.pathmover import IdentityPathMover
from openpathsampling.sample import Sample, SampleSet
from openpathsampling.shooting import UniformSelector
from openpathsampling.volume import CVDefinedVolume
import openpathsampling.engines.toy as toys
from .test_helpers import CallIdentity, raises_with_message_like
from .test_helpers import (assert_equal_array_array, items_equal,
                          make_1d_traj, assert_items_equal,
                          CalvinistDynamics, assert_same_items, A2BEnsemble)

#logging.getLogger('openpathsampling.pathmover').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.ensemble').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.netcdfplus').setLevel(logging.CRITICAL)


#logging.getLogger('openpathsampling.pathmover').propagate = False
#logging.getLogger('openpathsampling.initialization').propagate = False

class TestMakeListOfPairs(object):
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

def assert_sample_set_accepted(sample_set, results):
    for sample, result in zip(sample_set, results):
        assert_equal(sample.details.accepted, result)

def assert_subchanges_set_accepted(change, results):
    for ch, result in zip(change.subchanges, results):
        assert_equal(ch.accepted, result)

def assert_choice_of(result, choices):
    for choice in choices:
        if result is choice:
            return

    raise AssertionError("%s is not in list of choices [%s]" % (result, choices))

class TestPathMover(object):
    def setup(self):
        self.l1 = LengthEnsemble(1)
        self.l2 = LengthEnsemble(2)
        self.l3 = LengthEnsemble(3)
        self.s1 = Sample(replica=1, ensemble=self.l2)
        self.s2 = Sample(replica=2, ensemble=self.l1)
        self.s3 = Sample(replica=3, ensemble=self.l1)
        self.s4 = Sample(replica=2, ensemble=self.l3)
        self.sset = SampleSet([self.s1, self.s2, self.s3, self.s4])

    def test_legal_sample_set(self):
        assert_same_items(
            paths.PathMover.legal_sample_set(self.sset, ensembles=self.l1),
            [self.s2, self.s3]
        )
        assert_same_items(
            paths.PathMover.legal_sample_set(self.sset, ensembles=[self.l1]),
            [self.s2, self.s3]
        )

    def test_select_sample(self):
        for i in range(20):
            selected = PathMover.select_sample(self.sset)
            assert_choice_of(selected, [self.s1, self.s2, self.s3, self.s4])

    def test_is_ensemble_change_mover(self):
        pm = IdentityPathMover()
        assert_equal(pm.is_ensemble_change_mover, False)
        assert_equal(pm._is_ensemble_change_mover, None)
        pm._is_ensemble_change_mover = True
        assert_equal(pm.is_ensemble_change_mover, True)

    def test_is_canonical(self):
        pm = IdentityPathMover()
        assert_equal(pm.is_canonical, None)
        pm._is_canonical = True
        assert_equal(pm.is_canonical, True)


class TestShootingMover(object):
    def setup(self):
        self.dyn = CalvinistDynamics([-0.1, 0.1, 0.3, 0.5, 0.7,
                                      -0.1, 0.2, 0.4, 0.6, 0.8])
        op = FunctionCV("myid", f=lambda snap : snap.coordinates[0][0])
        self.stateA = CVDefinedVolume(op, -100, 0.0)
        self.stateB = CVDefinedVolume(op, 0.65, 100)
        self.tps = A2BEnsemble(self.stateA, self.stateB)
        init_traj = make_1d_traj(
            coordinates=[-0.1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7],
            velocities=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        )
        self.init_samp = SampleSet([Sample(
            trajectory=init_traj,
            replica=0,
            ensemble=self.tps
        )])

        integ = toys.LeapfrogVerletIntegrator(dt=0.1)
        pes = toys.LinearSlope(m=[0.0], c=[0.0])
        pes_AA = toys.LinearSlope(m=[10.0], c=[0.0])
        pes_BB = toys.LinearSlope(m=[-10.0], c=[0.0])
        topology = toys.Topology(n_spatial=1, masses=[1.0], pes=pes)
        self.toy_opts = {'integ': integ,
                         'n_frames_max': 1000,
                         'n_steps_per_frame': 1}
        # self.toy_engine: perfectly flat
        # self.toy_engine_AA: sloped to give AA trajectories
        # self.toy_engine_BB = sloped to give BB trajectories
        self.toy_engine = toys.Engine(options=self.toy_opts,
                                      topology=topology)
        self.toy_snap = toys.Snapshot(coordinates=np.array([[0.3]]),
                                      velocities=np.array([[0.1]]),
                                      engine=self.toy_engine)
        self.toy_traj = paths.Trajectory([
            toys.Snapshot(coordinates=np.array([[0.01*k - 0.005]]),
                          velocities=np.array([[0.1]]),
                          engine=self.toy_engine)
            for k in range(67)
        ])
        self.toy_samp = SampleSet([Sample(trajectory=self.toy_traj,
                                          replica=0,
                                          ensemble=self.tps)])
        # setup checks: keep around (commented) for debugging
        # assert_equal(self.stateA(self.toy_traj[0]), True)
        # assert_equal(self.stateA(self.toy_traj[1]), False)
        # assert_equal(self.stateB(self.toy_traj[-2]), False)
        # assert_equal(self.stateB(self.toy_traj[-1]), True)
        # assert_equal(self.tps(self.toy_traj), True)

class TestForwardShootMover(TestShootingMover):
    def test_move(self):
        mover = ForwardShootMover(
            ensemble=self.tps,
            selector=UniformSelector(),
            engine=self.dyn
        )
        self.dyn.initialized = True
        change = mover.move(self.init_samp)
        newsamp = self.init_samp.apply_samples(change)
        assert_equal(len(newsamp), 1)
        assert_equal(change.accepted, True)
        assert_equal(newsamp[0].ensemble(newsamp[0].trajectory), True)
        assert_equal(newsamp[0].trajectory, change.trials[0].trajectory)

    def test_move_toy_engine(self):
        mover = ForwardShootMover(
            ensemble=self.tps,
            selector=UniformSelector(),
            engine=self.toy_engine
        )
        change = mover.move(self.toy_samp)
        newsamp = self.toy_samp.apply_samples(change)
        assert_equal(len(newsamp), 1)
        assert_equal(change.accepted, True)
        newsamp.sanity_check()
        # the next two prove that we did a forward shot correctly
        assert_equal(self.toy_samp[0].trajectory[0],
                     newsamp[0].trajectory[0])
        assert_not_equal(self.toy_samp[0].trajectory[-1],
                         newsamp[0].trajectory[-1])


    def test_is_ensemble_change_mover(self):
        mover = ForwardShootMover(
            selector=UniformSelector(),
            ensemble=self.tps
        )

        assert_equal(mover.is_ensemble_change_mover, False)

class TestBackwardShootMover(TestShootingMover):
    def test_move(self):
        mover = BackwardShootMover(
            ensemble=self.tps,
            selector=UniformSelector(),
            engine=self.dyn
        )
        self.dyn.initialized = True
        change = mover.move(self.init_samp)
        newsamp = self.init_samp.apply_samples(change)
        assert_equal(len(newsamp), 1)
        assert_equal(len(self.init_samp[0].trajectory), 8)
        assert_not_equal(newsamp[0].trajectory[0],
                         self.init_samp[0].trajectory[0])
        assert_equal(newsamp[0].trajectory[-1],
                     self.init_samp[0].trajectory[-1])
        assert_true(len(newsamp[0].trajectory) <= 8)
        assert_equal(change.accepted, True)
        assert_equal(newsamp[0].ensemble(newsamp[0].trajectory), True)
        assert_equal(newsamp[0].trajectory, change.trials[0].trajectory)

    def test_move_toy_engine(self):
        mover = BackwardShootMover(
            ensemble=self.tps,
            selector=UniformSelector(),
            engine=self.toy_engine
        )
        change = mover.move(self.toy_samp)
        newsamp = self.toy_samp.apply_samples(change)
        assert_equal(len(newsamp), 1)
        assert_equal(change.accepted, True)
        newsamp.sanity_check()
        # the next two prove that we did a backward shot correctly
        assert_not_equal(self.toy_samp[0].trajectory[0],
                     newsamp[0].trajectory[0])
        assert_equal(self.toy_samp[0].trajectory[-1],
                         newsamp[0].trajectory[-1])


    def test_is_ensemble_change_mover(self):
        mover = BackwardShootMover(
            selector=UniformSelector(),
            ensemble=self.tps
        )
        assert_equal(mover.is_ensemble_change_mover, False)

class TestOneWayShootingMover(TestShootingMover):
    def test_mover_initialization(self):
        mover = OneWayShootingMover(
            ensemble=self.tps,
            selector=UniformSelector(),
            engine=self.dyn
        )
        assert_equal(len(mover.movers), 2)
        assert_equal(isinstance(mover, RandomChoiceMover), True)
        assert_equal(isinstance(mover, OneWayShootingMover), True)
        moverclasses = [m.__class__ for m in mover.movers]
        assert_equal(ForwardShootMover in moverclasses, True)
        assert_equal(BackwardShootMover in moverclasses, True)
        assert_equal(mover.engine, mover.movers[0].engine)
        assert_equal(mover.engine, mover.movers[1].engine)
        assert_equal(mover.selector, mover.movers[0].selector)
        assert_equal(mover.selector, mover.movers[1].selector)
        assert_equal(mover.ensemble, mover.movers[0].ensemble)
        assert_equal(mover.ensemble, mover.movers[1].ensemble)

class TestForwardFirstTwoWayShootingMover(TestShootingMover):
    _MoverType = ForwardFirstTwoWayShootingMover
    # this allows us to run the exact same tests for backward-first
    def test_run(self):
        mover = self._MoverType(
            ensemble=self.tps,
            selector=UniformSelector(),
            modifier=paths.NoModification(),
            engine=self.dyn
        )
        traj, details = mover._run(self.init_samp[0].trajectory, 4)
        assert_allclose(traj.xyz[:,0,0], [-0.1, 0.2, 0.4, 0.6, 0.8])
        assert_equal(list(details.keys()), ['modified_shooting_snapshot'])
        assert_equal(details['modified_shooting_snapshot'], traj[2])
        assert_not_in(details['modified_shooting_snapshot'],
                      self.init_samp[0].trajectory)

        traj, details = mover._run(self.init_samp[0].trajectory, 3)
        assert_allclose(traj.xyz[:,0,0], [-0.1, 0.1, 0.3, 0.5, 0.7])
        assert_equal(list(details.keys()), ['modified_shooting_snapshot'])
        assert_equal(details['modified_shooting_snapshot'], traj[2])
        assert_not_in(details['modified_shooting_snapshot'],
                      self.init_samp[0].trajectory)

    def test_run_toy(self):
        # mostly smoke test for toy engine integration
        mover = self._MoverType(
            ensemble=self.tps,
            selector=UniformSelector(),
            modifier=paths.NoModification(),
            engine=self.toy_engine
        )
        change = mover.move(self.toy_samp)
        assert_in(change.details.modified_shooting_snapshot,
                  change.trials[0].trajectory)
        assert_in(change.details.shooting_snapshot,
                  change.initial_trajectory)

    def _setup_early_reject(self, pair):
        zero_vel_traj = paths.Trajectory([
            toys.Snapshot(coordinates=np.array([[0.01*k - 0.005]]),
                          velocities=np.array([[0.0]]),
                          engine=self.toy_engine)
            for k in range(67)
        ])
        if pair == "AA":
            pes = toys.LinearSlope(m=[10.0], c=[0.0])
            traj = zero_vel_traj
        elif pair == "BB":
            pes = toys.LinearSlope(m=[-10.0], c=[0.0])
            traj = zero_vel_traj
        elif pair == "AB":
            pes = toys.LinearSlope(m=[0.0], c=[0.0])
            traj = self.init_samp[0].trajectory
        else:
            raise RuntimeError("unknown path type: " + pair)

        topology = toys.Topology(n_spatial=1, masses=[1.0], pes=pes)
        engine = toys.Engine(options=self.toy_opts, topology=topology)

        mapping = {'A': self.stateA, 'B': self.stateB}
        start = pair[0]
        end = pair[1]
        ensemble = paths.SequentialEnsemble([
            paths.AllInXEnsemble(mapping[start]) & paths.LengthEnsemble(1),
            paths.AllOutXEnsemble(self.stateA | self.stateB),
            paths.AllInXEnsemble(mapping[end]) & paths.LengthEnsemble(1)
        ])
        return (ensemble, engine, traj)

    def _test_early_reject(self, test_ensemble, path_types,
                           expected_rejections):
        for path_type in path_types:
            (ensemble, engine, traj) = self._setup_early_reject(path_type)
            # the ensemble returned above tells us what ensemble we expect
            # the trajectory to satisfy after the move, if *both* directions
            # were shot. (This is determined by the PES, which varies for
            # each path type.) We know whether we have early rejection based
            # on whether the trial trajectory satisfies that ensemble: if it
            # does, then we ran a full two-way shooting, and did not have
            # early rejection. If it does not, then we had early rejection.
            initial_sample_set = SampleSet([Sample(trajectory=traj,
                                                   replica=0,
                                                   ensemble=test_ensemble)])
            initial_sample_set.sanity_check()

            mover = self._MoverType(
                ensemble=test_ensemble,
                selector=UniformSelector(),
                modifier=paths.NoModification(),
                engine=engine
            )
            change = mover.move(initial_sample_set)

            expected_early_reject = path_type in expected_rejections
            ran_full_two_way = ensemble(change.trials[0].trajectory)
            assert_equal(expected_early_reject, not ran_full_two_way)

    def test_early_reject_tps(self):
        self._test_early_reject(test_ensemble=self.tps,
                                path_types=['AA', 'BB', 'AB'],
                                expected_rejections=['AA'])

    def test_early_reject_tis(self):
        both = self.stateA | self.stateB
        pseudo_tis_ensemble = paths.SequentialEnsemble([
            paths.AllInXEnsemble(self.stateA) & paths.LengthEnsemble(1),
            paths.AllOutXEnsemble(both),
            paths.AllInXEnsemble(both) & paths.LengthEnsemble(1)
        ])
        # when shooting forward-first, no TIS shot is ever rejected early
        self._test_early_reject(test_ensemble=pseudo_tis_ensemble,
                                path_types=['AA', 'BB', 'AB'],
                                expected_rejections=[])

    def test_sequential_shots(self):
        # make sure that, with no modification, the trajectory doesn't
        # change
        mover = self._MoverType(
            ensemble=self.tps,
            selector=UniformSelector(),
            modifier=paths.NoModification(),
            engine=self.toy_engine
        )
        change = mover.move(self.toy_samp)
        real_traj = change.trials[0].trajectory
        real_sample_set = SampleSet(change.trials)
        real_sample_set.sanity_check()

        new_change = mover.move(real_sample_set)
        new_traj = new_change.trials[0].trajectory
        assert_allclose(new_traj.xyz[:,0,0], real_traj.xyz[:,0,0])


class TestBackwardFirstTwoWayShootingMover(TestForwardFirstTwoWayShootingMover):
    _MoverType = BackwardFirstTwoWayShootingMover
    # runs the same tests as ForwardFirst
    def test_early_reject_tps(self):
        self._test_early_reject(test_ensemble=self.tps,
                                path_types=['AA', 'BB', 'AB'],
                                expected_rejections=['BB'])

    def test_early_reject_tis(self):
        both = self.stateA | self.stateB
        pseudo_tis_ensemble = paths.SequentialEnsemble([
            paths.AllInXEnsemble(self.stateA) & paths.LengthEnsemble(1),
            paths.AllOutXEnsemble(both),
            paths.AllInXEnsemble(both) & paths.LengthEnsemble(1)
        ])
        # when shooting forward-first, no TIS shot is ever rejected early
        self._test_early_reject(test_ensemble=pseudo_tis_ensemble,
                                path_types=['AA', 'BB', 'AB'],
                                expected_rejections=['BB'])


class TestTwoWayShootingMover(TestShootingMover):
    def test_properties(self):
        selector = UniformSelector()
        modifier = paths.NoModification()
        mover = TwoWayShootingMover(
            ensemble=self.tps,
            selector=selector,
            modifier=modifier,
            engine=self.dyn
        )
        assert_equal(mover.ensemble, self.tps)
        assert_equal(mover.selector, selector)
        assert_equal(mover.modifier, modifier)

    def test_to_dict_from_dict(self):
        mover = TwoWayShootingMover(
            ensemble=self.tps,
            selector=UniformSelector(),
            modifier=paths.NoModification(),
            engine=self.dyn
        )
        dct = mover.to_dict()
        new_mover = mover.from_dict(dct)
        assert_equal(mover.movers, new_mover.movers)
        assert_equal(mover.ensemble, new_mover.ensemble)
        assert_equal(mover.selector, new_mover.selector)
        assert_equal(mover.modifier, new_mover.modifier)


class TestPathReversalMover(object):
    def setup(self):
        op = FunctionCV("myid", f=lambda snap :
                             snap.coordinates[0][0])

        volA = CVDefinedVolume(op, -100, 0.0)
        volB = CVDefinedVolume(op, 1.0, 100)
        volX = CVDefinedVolume(op, -100, 0.25)
        self.tis = paths.TISEnsemble(volA, volB, volX, orderparameter=op,
                                     lambda_i=0.25)
        self.move = PathReversalMover(ensemble=self.tis)
        self.op = op

    def test_is_ensemble_change_mover(self):
        assert_equal(self.move.is_ensemble_change_mover, False)

    def test_AXA_path(self):
        trajAXA = make_1d_traj(coordinates=[-0.1, 0.75, -0.6],
                               velocities=[0.1, 0.05, -0.05])
        assert_equal(self.tis(trajAXA), True)
        sampAXA = Sample(trajectory=trajAXA,
                         ensemble=self.tis,
                         replica=0)
        gs_AXA = SampleSet([sampAXA])
        change = self.move.move(gs_AXA)
        assert_equal(change.accepted, True)

    def test_A_A_path(self):
        trajA_A = make_1d_traj(coordinates=[-0.3, 0.1, -0.4])
        sampA_A = Sample(trajectory=trajA_A,
                         ensemble=self.tis,
                         replica=0)
        gs_A_A = SampleSet([sampA_A])
        change = self.move.move(gs_A_A)
        assert_equal(change.accepted, False)


    def test_AB_path(self):
        trajAXB = make_1d_traj(coordinates=[-0.2, 0.75, 1.8])
        sampAXB = Sample(trajectory=trajAXB,
                         ensemble=self.tis,
                         replica=0)
        gs_AXB = SampleSet([sampAXB])
        change = self.move.move(gs_AXB)
        assert_equal(change.accepted, False)


    def test_BA_path(self):
        trajBXA = make_1d_traj(coordinates=[1.2, 0.7, -0.25])
        sampBXA = Sample(trajectory=trajBXA,
                         ensemble=self.tis,
                         replica=0)
        gs_BXA = SampleSet([sampBXA])
        change = self.move.move(gs_BXA)
        # print [[v.coordinates[0] for v in t.trajectory] for t in change.trials]
        assert_equal(change.accepted, True)


class TestReplicaIDChangeMover(object):
    def setup(self):
        pass

    def test_replica_in_sample_set(self):
        raise SkipTest

    def test_replica_not_in_sample_set(self):
        raise SkipTest


class TestReplicaExchangeMover(object):
    def setup(self):
        op = FunctionCV("myid", f=lambda snap :
                             snap.coordinates[0][0])

        state1 = CVDefinedVolume(op, -100, 0.0)
        state2 = CVDefinedVolume(op, 1, 100)
        volA = CVDefinedVolume(op, -100, 0.25)
        volB = CVDefinedVolume(op, -100, 0.50)
        self.tisA = paths.TISEnsemble(state1, state2, volA,
                                      orderparameter=op, lambda_i=0.25)
        self.tisB = paths.TISEnsemble(state1, state2, volB,
                                      orderparameter=op, lambda_i=0.50)
        self.traj0 = make_1d_traj([-0.1, 0.2, 0.3, 0.1, -0.2])
        self.traj1 = make_1d_traj([-0.1, 0.1, 0.4, 0.6, 0.3, 0.2, -0.15])
        self.traj2 = make_1d_traj([-0.1, 0.2, 0.3, 0.7, 0.6, 0.4, 0.1, -0.15])
        self.sampA0 = Sample(replica=0, trajectory=self.traj0, ensemble=self.tisA)
        self.sampB1 = Sample(replica=1, trajectory=self.traj1, ensemble=self.tisB)
        self.sampA2 = Sample(replica=2, trajectory=self.traj2, ensemble=self.tisA)
        self.gs_B1A2 = SampleSet([self.sampB1, self.sampA2])
        self.gs_A0B1 = SampleSet([self.sampA0, self.sampB1])

    def test_is_ensemble_change_mover(self):
        repex_AB = ReplicaExchangeMover(
            ensemble1=self.tisA,
            ensemble2=self.tisB
        )
        assert_equal(repex_AB.is_ensemble_change_mover, True)

    def test_repex_ens_rej(self):
        repex_AB = ReplicaExchangeMover(
            ensemble1=self.tisA,
            ensemble2=self.tisB
        )
        old_sset = self.gs_A0B1
        repex_change = repex_AB.move(old_sset)
        samples = repex_change.results

        assert_equal(len(repex_change.results), 0) # since rejected

        samples_A0B1_ens = repex_change.trials
        assert_equal(len(samples_A0B1_ens), 2)
        assert_equal(repex_change.accepted, False)

        new_sset = old_sset.apply_samples(samples)

        assert_equal(new_sset[0].trajectory, old_sset[0].trajectory)
        assert_equal(new_sset[1].trajectory, old_sset[1].trajectory)

        A0 = [s for s in samples_A0B1_ens if s.ensemble==self.tisA]
        assert_equal(len(A0), 1)
        assert_equal(A0[0].trajectory, self.traj1)
        B1 = [s for s in samples_A0B1_ens if s.ensemble==self.tisB]
        assert_equal(len(B1), 1)
        assert_equal(B1[0].trajectory, self.traj0)

    def test_repex_ens_acc(self):
        repex_12 = ReplicaExchangeMover(
            ensemble1=self.tisA,
            ensemble2=self.tisB
        )
        old_sset = self.gs_B1A2
        samples_B2A1_rep = repex_12.move(old_sset)
        change = samples_B2A1_rep
        samples = change.results
        assert_equal(len(samples), 2)

        assert_equal(change.accepted, True)

        new_sset = old_sset.apply_samples(samples)

        assert_equal(new_sset[1].trajectory, old_sset[1].trajectory)
        assert_equal(new_sset[2].trajectory, old_sset[2].trajectory)

        B2 = [s for s in samples if s.ensemble==self.tisB]
        assert_equal(len(B2), 1)
        assert_equal(B2[0].trajectory, self.traj2)
        A1 = [s for s in samples if s.ensemble==self.tisA]
        assert_equal(len(A1), 1)
        assert_equal(A1[0].trajectory, self.traj1)


class TestRandomChoiceMover(object):
    def setup(self):
        traj = Trajectory([-0.5, 0.7, 1.1])
        op = CallIdentity()
        volA = CVDefinedVolume(op, -100, 0.0)
        volB = CVDefinedVolume(op, 1.0, 100)
        volX = CVDefinedVolume(op, -100, 0.25)
        self.tis = paths.TISEnsemble(volA, volB, volX)
        self.tps = A2BEnsemble(volA, volB)
        self.len3 = LengthEnsemble(3)
        self.init_samp = SampleSet([Sample(trajectory=traj,
                                           ensemble=self.len3,
                                           replica=0)])
        self.hop_to_tis = EnsembleHopMover(
            ensemble=self.len3,
            target_ensemble=self.tis
        )
        self.hop_to_tps = EnsembleHopMover(
            ensemble=self.len3,
            target_ensemble=self.tps
        )
        self.mover = RandomChoiceMover([self.hop_to_tis, self.hop_to_tps])

    def test_is_ensemble_change_mover(self):
        assert_equal(self.mover.is_ensemble_change_mover, True)

    def test_is_canonical(self):
        for t in range(20):
            change = self.mover.move(self.init_samp)
            assert_not_equal(change.canonical.mover, self.mover)
            canonical_submovers = 0
            for submover in self.mover.movers:
                if change.canonical.mover is submover:
                    canonical_submovers += 1
            assert_equal(canonical_submovers, 1)

    def test_random_choice(self):
        # test that both get selected, but that we always return only one
        # sample
        count = {}
        for t in range(100):
            change = self.mover.move(self.init_samp)
            assert_equal(len(change.results), 1)
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


class TestRandomAllowedChoiceMover(object):
    def setup(self):
        self.dyn = CalvinistDynamics([-0.1, 0.1, 0.3, 0.5, 0.7,
                                      -0.1, 0.2, 0.4, 0.6, 0.8,
                                     ])
        self.dyn.initialized = True
        # SampleMover.engine = self.dyn
        op = FunctionCV("myid", f=lambda snap :
                             snap.coordinates[0][0])
        stateA = CVDefinedVolume(op, -100, 0.0)
        stateB = CVDefinedVolume(op, 0.65, 100)
        volX = CVDefinedVolume(op, -100, 0.25)
        volY = CVDefinedVolume(op, -100, 0.40)
        self.ens1 = paths.TISEnsemble(stateA, stateB, volX, op)
        self.ens2 = paths.TISEnsemble(stateA, stateB, volY, op)
        init_traj1 = make_1d_traj(
            coordinates=[-0.1, 0.1, 0.2, 0.3, 0.24, 0.15, 0.06, -0.07],
            velocities=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        )
        init_traj2 = make_1d_traj(
            coordinates=[-0.1, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7],
            velocities=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
        )
        self.samp1 = Sample(trajectory=init_traj1, replica=0,
                            ensemble=self.ens1)
        self.samp2 = Sample(trajectory=init_traj2, replica=1,
                            ensemble=self.ens2)

        self.shooter = ForwardShootMover(selector=UniformSelector(),
                                         ensemble=self.ens2,
                                         engine=self.dyn)
        self.pathrev = PathReversalMover(ensemble=self.ens1)

        ens_dict = {self.ens1 : self.pathrev, self.ens2 : self.shooter}
        # self.mover = EnsembleDictionaryMover(ens_dict)
        self.mover = RandomAllowedChoiceMover([self.shooter, self.pathrev])

    def test_move_single_replica(self):
        sample_set = SampleSet([self.samp1])
        change = self.mover.move(sample_set)
        subchange = change.subchange
        assert_equal(subchange.mover, self.pathrev)
        assert_equal(subchange.accepted, True)
        assert_equal(change.accepted, True)
        assert_equal(len(subchange.samples), 1)

        sample_set = SampleSet([self.samp2])
        change = self.mover.move(sample_set)
        subchange = change.subchange
        assert_equal(subchange.mover, self.shooter)
        assert_equal(subchange.accepted, True)
        assert_equal(change.accepted, True)


    def test_move_multiple_replicas(self):
        sample_set = SampleSet([self.samp1, self.samp2])
        count = {}
        for i in range(100):
            change = self.mover.move(sample_set)
            subchange = change.subchange
            assert_equal(change.accepted, True)
            assert_equal(subchange.accepted, True)
            assert_equal(len(subchange.samples), 1)
            ens = subchange.trials[0].ensemble
            try:
                count[ens] += 1
            except KeyError:
                count[ens] = 1
            if ens == self.ens1:
                assert_equal(subchange.mover, self.pathrev)
            elif ens == self.ens2:
                assert_equal(subchange.mover, self.shooter)
            else:
                raise AssertionError("Resulting mover unknown!")
        assert_equal(set(count.keys()), set([self.ens1, self.ens2]))


    def test_move_multiple_replicas_weighted_ensembles(self):
        sample_set = SampleSet([self.samp1, self.samp2])
        ens_dict = {self.ens1 : self.pathrev, self.ens2 : self.shooter}
        # weighted_mover = EnsembleDictionaryMover(ens_dict, [1.0, 2.0])
        weighted_mover = RandomAllowedChoiceMover([self.pathrev,
                                                   self.shooter], [1.0, 2.0])
        count = {}
        for i in range(100):
            change = weighted_mover.move(sample_set)
            subchange = change.subchange
            assert_equal(change.accepted, True)
            assert_equal(subchange.accepted, True)
            assert_equal(len(subchange.samples), 1)
            ens = subchange.trials[0].ensemble
            try:
                count[ens] += 1
            except KeyError:
                count[ens] = 1
            if ens == self.ens1:
                assert_equal(subchange.mover, self.pathrev)
            elif ens == self.ens2:
                assert_equal(subchange.mover, self.shooter)
            else:
                raise AssertionError("Resulting mover unknown!")
        assert_equal(set(count.keys()), set([self.ens1, self.ens2]))
        try:
            assert(count[self.ens1] < count[self.ens2])
        except AssertionError:
            raise AssertionError("Not true: "+str(count[self.ens1]) + " < "
                                 + str(count[self.ens2]))


class TestSequentialMover(object):
    def setup(self):
        traj = Trajectory([-0.5, 0.7, 1.1])
        op = CallIdentity()
        volA = CVDefinedVolume(op, -100, 0.0)
        volB = CVDefinedVolume(op, 1.0, 100)
        volX = CVDefinedVolume(op, -100, 0.25)
        tis = paths.TISEnsemble(volA, volB, volX)
        tps = A2BEnsemble(volA, volB)
        len3 = LengthEnsemble(3)
        len2 = LengthEnsemble(2)
        self.hop_to_tis = RandomAllowedChoiceMover(
            [EnsembleHopMover(*ens) for ens in [[tis, tis],
                 [tps, tis],
                 [len3, tis],
                 [len2, tis]]]
        )
        self.hop_to_tps = RandomAllowedChoiceMover(
            [EnsembleHopMover(*ens) for ens in [[tis, tps],
                 [tps, tps],
                 [len3, tps],
                 [len2, tps]]]
        )
        self.hop_to_len3 = RandomAllowedChoiceMover(
            [EnsembleHopMover(*ens) for ens in [[tis, len3],
                 [tps, len3],
                 [len3, len3],
                 [len2, len3]]]
        )
        self.hop_to_len2 = RandomAllowedChoiceMover(
            [EnsembleHopMover(*ens) for ens in [[tis, len2],
                 [tps, len2],
                 [len3, len2],
                 [len2, len2]]]
        )
        self.init_sample = Sample(trajectory=traj,
                                  ensemble=len3,
                                  replica=0)
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

    def test_is_ensemble_change_mover(self):
        move = SequentialMover(movers=self.everything_accepted_movers)
        assert_equal(move.is_ensemble_change_mover, True)

    def test_everything_accepted(self):
        move = SequentialMover(movers=self.everything_accepted_movers)
        gs = SampleSet(self.init_sample)
        change = move.move(gs)
        samples = change.results
        assert_equal(len(samples), 3)
        for subchange in change:
            assert_equal(subchange.accepted, True)
        gs = gs.apply_samples(change)
        assert_equal(gs[0].ensemble, self.tps)

    def test_first_rejected(self):
        move = SequentialMover(movers=self.first_rejected_movers)
        gs = SampleSet(self.init_sample)
        change = move.move(gs)
        samples = change.results
        # @DWHS: This should have two samples since two are accepted
        # and thus applied
        assert_equal(len(samples), 2)
        assert_equal(change[0].accepted, False)
        assert_equal(change[1].accepted, True)
        assert_equal(change[2].accepted, True)
        gs = gs.apply_samples(change)
        assert_equal(gs[0].ensemble, self.tps)

    def test_last_rejected(self):
        move = SequentialMover(movers=self.last_rejected_movers)
        gs = SampleSet(self.init_sample)
        change = move.move(gs)
        samples = change.results
        assert_equal(len(samples), 2)
        # @DWHS: I think if the last is rejected then there should only be two
        # samples to be used, since the last one is not accepted and thus
        # discarded (does not mean that it is not stored!!!)
        assert_equal(change[0].accepted, True)
        assert_equal(change[1].accepted, True)
        assert_equal(change[2].accepted, False)
        gs = gs.apply_samples(change)
        assert_equal(gs[0].ensemble, self.tps)

    def test_restricted_by_replica(self):
        raise SkipTest

    def test_restricted_by_ensemble(self):
        raise SkipTest

class TestPartialAcceptanceSequentialMover(TestSequentialMover):
    def test_everything_accepted(self):
        move = PartialAcceptanceSequentialMover(movers=self.everything_accepted_movers)
        gs = SampleSet(self.init_sample)
        change = move.move(gs)
        samples = change.results
        assert_equal(len(samples), 3)
        for subchange in change:
            assert_equal(subchange.accepted, True)
        assert_equal(len(change.trials,),3)
        gs = gs.apply_samples(change)
        assert_equal(gs[0].ensemble, self.tps)

    def test_first_rejected(self):
        move = PartialAcceptanceSequentialMover(movers=self.first_rejected_movers)
        gs = SampleSet(self.init_sample)
        change = move.move(gs)
        samples = change.results
        # returns zero sample since even the first is rejected
        # the first one is still stored
        assert_equal(len(samples), 0)
        allsamp = change.trials
        assert_equal(len(allsamp), 1)
        assert_equal(change[0].accepted, False)
        gs = gs.apply_samples(change)
        assert_equal(gs[0].ensemble, self.len3)

    def test_last_rejected(self):
        move = PartialAcceptanceSequentialMover(movers=self.last_rejected_movers)
        gs = SampleSet(self.init_sample)
        change = move.move(gs)
        samples = change.results
        # @see above, this should return 2 samples. Important the third is
        # still run!
        assert_equal(len(samples), 2)
        allsamp = change.trials
        assert_equal(len(allsamp), 3)

        assert_equal(change[0].accepted, True)
        assert_equal(change[1].accepted, True)
        assert_equal(change[2].accepted, False)
        gs = gs.apply_samples(change)
        assert_equal(gs[0].ensemble, self.tps)

    def test_restricted_by_replica(self):
        raise SkipTest

    def test_restricted_by_ensemble(self):
        raise SkipTest

class TestConditionalSequentialMover(TestSequentialMover):
    def test_everything_accepted(self):
        move = ConditionalSequentialMover(movers=self.everything_accepted_movers)
        gs = SampleSet(self.init_sample)
        change = move.move(gs)
        samples = change.results
        assert_equal(len(samples), 3)
        for ch in change:
            assert_equal(change.accepted, True)
        gs = gs.apply_samples(change)
        assert_equal(gs[0].ensemble, self.tps)

    def test_first_rejected(self):
        move = ConditionalSequentialMover(movers=self.first_rejected_movers)
        gs = SampleSet(self.init_sample)
        change = move.move(gs)
        samples = change.results
        # should be zero since the move is completely rejected
        assert_equal(len(samples), 0)
        allsamp = change.trials
        assert_equal(len(allsamp), 1)
        assert_equal(change[0].accepted, False)
        gs = gs.apply_samples(change)
        assert_equal(gs[0].ensemble, self.len3)

    def test_last_rejected(self):
        move = ConditionalSequentialMover(movers=self.last_rejected_movers)
        gs = SampleSet(self.init_sample)
        change = move.move(gs)
        samples = change.results
        # number of accepted samples is 0 for this type of mover
        assert_equal(len(samples), 0)
        allsamp = change.trials
        assert_equal(len(allsamp), 3)

        # check here if last actual samples was false
        # this actually allows to see later if the single samples were
        # accepted or not, even from the change without loading samples
        assert_equal(change[0].accepted, True)
        assert_equal(change[1].accepted, True)
        assert_equal(change[2].accepted, False)
        gs = gs.apply_samples(change)
        assert_equal(gs[0].ensemble, self.len3)

    def test_restricted_by_replica(self):
        raise SkipTest

    def test_restricted_by_ensemble(self):
        raise SkipTest

class SubtrajectorySelectTester(object):

    def setup(self):
        op = CallIdentity()
        vol = paths.CVDefinedVolume(op, -0.5, 0.5)
        inX = paths.AllInXEnsemble(vol)
        outX = paths.AllOutXEnsemble(vol)
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
            ensemble=self.ensemble,
            trajectory=self.traj_with_3_subtrajs
        ))

    def test_paths_in_ensemble(self):
        # more a test of SequentialEnsemble, but also a test of sanity
        # before the real tests
        assert_equal(self.ensemble(self.traj_with_3_subtrajs), True)
        assert_equal(self.subensemble(self.subtraj0), True)
        assert_equal(self.subensemble(self.subtraj1), True)
        assert_equal(self.subensemble(self.subtraj2), True)

class TestRandomSubtrajectorySelectMover(SubtrajectorySelectTester):
    def test_accepts_all(self):
        mover = RandomSubtrajectorySelectMover(
            ensemble=self.ensemble,
            sub_ensemble=self.subensemble
        )
        found = {}
        for t in range(100):
            change = mover.move(self.gs)
            samples = change.results
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

    def test_is_ensemble_change_mover(self):
        mover = RandomSubtrajectorySelectMover(ensemble=1, sub_ensemble=1)
        assert_equal(mover.is_ensemble_change_mover, True)


    def test_nl_fails(self):
        raise SkipTest

    def test_nothing_allowed(self):
        mover = RandomSubtrajectorySelectMover(
            ensemble=self.ensemble,
            sub_ensemble=self.subensemble
        )
        traj_with_no_subtrajs = Trajectory([0.0, 0.0, 0.0])
        self.gs[0].trajectory = traj_with_no_subtrajs
        change = mover.move(self.gs)
        samples = change.results
        assert_equal(len(samples), 0)
        # print change.samples
        assert_equal(len(change.samples), 0)

class TestFirstSubtrajectorySelectMover(SubtrajectorySelectTester):
    def test_move(self):
        mover = FirstSubtrajectorySelectMover(
            ensemble=self.ensemble,
            sub_ensemble=self.subensemble
        )
        change = mover.move(self.gs)
        samples = change.results
        assert_equal(len(samples), 1)
        assert_equal(self.subensemble, samples[0].ensemble)
        assert_equal(self.subensemble(samples[0].trajectory), True)
        assert_equal(self.ensemble(samples[0].trajectory), False)
        assert_equal(samples[0].trajectory, self.subtraj0)

class TestFinalSubtrajectorySelectMover(SubtrajectorySelectTester):
    def test_move(self):
        mover = FinalSubtrajectorySelectMover(
            ensemble=self.ensemble,
            sub_ensemble=self.subensemble
        )
        change = mover.move(self.gs)
        samples = change.results
        assert_equal(len(samples), 1)
        assert_equal(self.subensemble, samples[0].ensemble)
        assert_equal(self.subensemble(samples[0].trajectory), True)
        assert_equal(self.ensemble(samples[0].trajectory), False)
        assert_equal(samples[0].trajectory, self.subtraj2)

# class TestForceEnsembleChangeMover(object):
#     def setup(self):
#         traj = Trajectory([-0.5, 0.7, 1.1])
#         op = CallIdentity()
#         volA = CVDefinedVolume(op, -100, 0.0)
#         volB = CVDefinedVolume(op, 1.0, 100)
#         volX = CVDefinedVolume(op, -100, 0.25)
#         self.tis = paths.TISEnsemble(volA, volB, volX)
#         self.len3 = LengthEnsemble(3)
#         self.len2 = LengthEnsemble(2)
#         self.gs = SampleSet(Sample(
#             trajectory=traj,
#             ensemble=self.tis,
#             replica=0
#         ))
#
#     def test_in_ensemble(self):
#         mover = ForceEnsembleChangeMover(ensembles=[[self.tis, self.len3]])
#         change = mover.move(self.gs)
#         samples = change.results
#         assert_equal(change.details.initial_ensemble(samples[0].trajectory),
#                      True)
#         assert_equal(samples[0].ensemble(samples[0].trajectory), True)
#         assert_equal(samples[0].ensemble, self.len3)
#
#     def test_not_in_ensemble(self):
#         mover = ForceEnsembleChangeMover(ensembles=[[self.tis, self.len2]])
#         change = mover.move(self.gs)
#         samples = change.results
#         assert_equal(change.details.initial_ensemble(samples[0].trajectory),
#                      True)
#         assert_equal(samples[0].ensemble, self.len2)
#         assert_equal(samples[0].ensemble(samples[0].trajectory), False)

class TestMinusMover(object):
    def setup(self):
        op = FunctionCV("myid", f=lambda snap :
                             snap.coordinates[0][0])

        volA = CVDefinedVolume(op, -100, 0.0)
        volB = CVDefinedVolume(op, 1.0, 100)
        volX = CVDefinedVolume(op, -100, 0.25)
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
        # SampleMover.engine = self.dyn
        self.dyn.initialized = True
        self.innermost = paths.TISEnsemble(volA, volB, volX)
        self.minus = paths.MinusInterfaceEnsemble(volA, volX)
        self.mover = MinusMover(
            minus_ensemble=self.minus,
            innermost_ensembles=self.innermost,
            engine=self.dyn
        )
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

    def test_is_ensemble_change_mover(self):
        assert_equal(self.mover.is_ensemble_change_mover, True)

    def test_is_canonical(self):
        assert_equal(self.mover.is_canonical, True)

    def test_setup_sanity(self):
        # sanity checks to make sure that what we set up makes sense
        assert_equal(self.minus_sample.ensemble(self.minus_sample.trajectory),
                    True)
        first_subtraj = FirstSubtrajectorySelectMover(
            ensemble=self.minus,
            sub_ensemble=self.minus._segment_ensemble
        )
        change = first_subtraj.move(SampleSet(self.minus_sample))
        samples = change.results
        assert_equal(samples[0].ensemble(samples[0].trajectory), True)
        final_subtraj = FinalSubtrajectorySelectMover(
            ensemble=self.minus,
            sub_ensemble=self.minus._segment_ensemble
        )
        change = final_subtraj.move(SampleSet(self.minus_sample))
        samples = change.results
        assert_equal(samples[0].ensemble(samples[0].trajectory), True)
        assert_equal(samples[0].ensemble, self.minus._segment_ensemble)
        assert_equal(self.mover.engine, self.dyn)


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
            change = self.mover.move(gs)
            assert change.details.segment_swap_samples is not None
            assert change.details.extension_trajectory is not None
            samples = change.results
            sub_samples = change.subchange.subchange.results
            assert_equal(len(samples), 2)
            assert_equal(len(sub_samples), 4)
            s_inner = [s for s in sub_samples if s.ensemble==self.innermost]
            s_minus = [s for s in sub_samples if s.ensemble==self.minus]
            s_sub = [s for s in sub_samples if s.ensemble==self.minus._segment_ensemble]
            assert_equal(len(s_inner), 1)
            assert_equal(len(s_minus), 1)
            assert_equal(len(s_sub), 2)

            for c in change:
                assert_equal(c.accepted, True)

            assert_equal(change.canonical.mover, self.mover)

            key = ""
            s_inner0_xvals = [s.coordinates[0,0] for s in s_inner[0].trajectory]
            if items_equal(s_inner0_xvals, self.first_segment):
                key += "1"
            elif items_equal(s_inner0_xvals, self.second_segment):
                key += "2"
            else:
                print("s_inner0_xvals:", s_inner0_xvals)
                raise RuntimeError("Chosen segment neither first nor last!")

            # final sample s_minus is accepted
            s_minus_xvals = [s.coordinates[0,0] for s in s_minus[-1].trajectory]
            if items_equal(s_minus_xvals, extend_forward):
                key += "f"
            elif items_equal(s_minus_xvals, extend_backward):
                key += "b"
            else:
                print("s_minus_xvals:", s_minus_xvals)
                raise RuntimeError("Unexpected minus extension result!")

            try:
                seg_dir[key] += 1
            except KeyError:
                seg_dir[key] = 1
        assert_equal(len(list(seg_dir.keys())), 4)

    def test_repex_fails_other_ensemble(self):
        innermost_other_ensemble = make_1d_traj([-0.11, 0.1, -0.12])
        samp_other_ensemble = Sample(
            replica=0,
            trajectory=innermost_other_ensemble,
            ensemble=self.innermost
        )
        gs = SampleSet([samp_other_ensemble, self.minus_sample])

        change = self.mover.move(gs)
        assert_equal(len(change.trials), 1)
        assert change.details.segment_swap_samples is not None
        assert change.details.extension_trajectory is None

        sub = change.subchange.subchange
        assert_equal(self.innermost(innermost_other_ensemble), False)
        assert_equal(sub[0].accepted, True)
        assert_equal(sub[1].accepted, False)
        assert_equal(len(sub.trials), 3) # stop after failed repex
        # only one sample which is not a segment

    def test_repex_fails_innermost_crosses_state(self):
        innermost_crosses_to_state = make_1d_traj([-0.11, 0.5, 1.8])
        samp_crosses_to_state = Sample(
            replica=0,
            trajectory=innermost_crosses_to_state,
            ensemble=self.innermost
        )
        gs = SampleSet([samp_crosses_to_state, self.minus_sample])

        change = self.mover.move(gs)
        assert_equal(len(change.trials), 1) # stop after failed repex
        assert change.details.segment_swap_samples is not None
        assert change.details.extension_trajectory is None

        sub = change.subchange.subchange
        assert_equal(self.innermost(innermost_crosses_to_state), True)
        assert_equal(len(sub.trials), 3) # stop after failed repex
        assert_subchanges_set_accepted(sub, [True, False, False])

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

        change = self.mover.move(gs)
        assert change.details.segment_swap_samples is not None
        assert change.details.extension_trajectory is None
        sub = change.subchange.subchange
        assert_equal(len(sub.trials), 3)  # stop after failed repex
        assert_equal(len(change.trials), 1)
        assert_subchanges_set_accepted(sub, [True, False, False])

    def test_extension_fails(self):

        # we use `stop` and not `fail` for max length
        self.dyn.options['on_max_length'] = 'stop'

        innermost_bad_extension = [-0.25, 0.1, 0.5, 0.1, -0.25]
        traj_bad_extension = make_1d_traj(innermost_bad_extension, [1.0]*5)
        samp_bad_extension = Sample(
            replica=0,
            trajectory=traj_bad_extension,
            ensemble=self.innermost
        )

        assert_equal(self.innermost(traj_bad_extension), True)

        gs = SampleSet([self.minus_sample, samp_bad_extension])
        change = self.mover.move(gs)
        assert_equal(change.accepted, False) # whole minus has failed

        sub = change.subchange.subchange
        assert_equal(len(sub.trials), 4)
        assert change.details.segment_swap_samples is not None
        assert change.details.extension_trajectory is not None

        # after filtering there are only 2 trials
        assert_equal(len(change.trials), 2)

        assert_subchanges_set_accepted(sub, [True] * 2 + [False])
        # first two work and the extension fails
        # this only happens due to length

        assert_equal(
            len(sub[-1][0].trials[0].trajectory),
            len(traj_bad_extension)+self.dyn.n_frames_max-1
        )

        self.dyn.options['on_max_length'] = 'fail'


class TestSingleReplicaMinusMover(object):
    def setup(self):
        op = FunctionCV("myid", f=lambda snap :
                             snap.coordinates[0][0])

        volA = CVDefinedVolume(op, -100, 0.0)
        volB = CVDefinedVolume(op, 1.0, 100)
        volX = CVDefinedVolume(op, -100, 0.25)
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
        # SampleMover.engine = self.dyn
        self.dyn.initialized = True
        self.innermost = paths.TISEnsemble(volA, volB, volX)
        self.minus = paths.MinusInterfaceEnsemble(volA, volX)
        self.mover = SingleReplicaMinusMover(
            minus_ensemble=self.minus,
            innermost_ensembles=self.innermost,
            engine=self.dyn
        )
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

    def test_is_ensemble_change_mover(self):
        assert_equal(self.mover.is_ensemble_change_mover, True)

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

        output_forward = [-0.12, 0.12, 0.32, -0.131]
        output_backward = [-0.13, 0.13, 0.33, -0.11]

        seg_dir = {}
        for i in range(100):
            change = self.mover.move(gs)
            samples = change.results
            sub_samples = change.subchange.subchange.results
            assert_equal(len(samples), 1)
            assert_equal(len(sub_samples), 4)
            s_inner = [s for s in sub_samples if s.ensemble==self.innermost]
            s_minus = [s for s in sub_samples if s.ensemble==self.minus]
            s_seg = [s for s in sub_samples if s.ensemble==self.minus._segment_ensemble]
            assert_equal(len(s_inner), 1) # this is the output
            assert_equal(len(s_minus), 1) # this is the minus version
            assert_equal(len(s_seg), 2) # first the selected, then the final

            for c in change:
               assert_equal(c.accepted, True)

            assert_equal(change.canonical.mover, self.mover)

            key = ""
            s_seg0_xvals = [s.coordinates[0,0] for s in s_seg[0].trajectory]
            if items_equal(s_seg0_xvals, self.list_innermost):
                key += "0"
            else:
                print("s_seg0_xvals:", s_seg0_xvals)
                raise RuntimeError("Chosen segment neither first nor last!")

            # s_minus is the intermediate
            s_minus_xvals = [s.coordinates[0,0] for s in s_minus[-1].trajectory]
            if items_equal(s_minus_xvals, extend_forward):
                key += "f"
                assert_equal(
                    [s.coordinates[0,0] for s in s_inner[0].trajectory],
                    output_forward
                )
            elif items_equal(s_minus_xvals, extend_backward):
                key += "b"
                assert_equal(
                    [s.coordinates[0,0] for s in s_inner[0].trajectory],
                    output_backward
                )
            else:
                print("s_minus_xvals:", s_minus_xvals)
                raise RuntimeError("Unexpected minus extension result!")


            try:
                seg_dir[key] += 1
            except KeyError:
                seg_dir[key] = 1
        assert_equal(len(list(seg_dir.keys())), 2)

    def test_first_hop_fails(self):
        crossing_traj = make_1d_traj([-0.11, 0.11, 0.31, 1.01], [1.0]*4)
        crossing_samp = Sample(replica=0, trajectory=crossing_traj,
                               ensemble=self.innermost)
        gs = SampleSet([crossing_samp])
        gs.sanity_check()

        change = self.mover.move(gs)
        assert_equal(change.accepted, False)
        assert_equal(len(change.results), 0)
        sub_trials = change.subchange.subchange.subchange.trials
        assert_equal(len(sub_trials), 1)
        assert_equal(sub_trials[0].trajectory, crossing_traj)
        assert_equal(sub_trials[0].ensemble, self.minus._segment_ensemble)

    def test_extension_fails(self):
        # we use `stop` and not `fail` for max length
        self.dyn.options['on_max_length'] = 'stop'

        innermost_bad_extension = [-0.25, 0.1, 0.5, 0.1, -0.25]
        traj_bad_extension = make_1d_traj(innermost_bad_extension, [1.0]*5)
        samp_bad_extension = Sample(
            replica=0,
            trajectory=traj_bad_extension,
            ensemble=self.innermost
        )

        assert_equal(self.innermost(traj_bad_extension), True)

        gs = SampleSet([self.minus_sample, samp_bad_extension])
        change = self.mover.move(gs)
        assert_equal(change.accepted, False) # whole minus has failed

        #     Minus : Filter  :ChooseFB : CondSeq
        sub = change.subchange.subchange.subchange
        assert_equal(len(sub.trials), 2)
        assert_equal(len(change.trials), 0) # no trials survive filtering
        assert_subchanges_set_accepted(sub, [True, False])

        # first two work and the extension fails
        # this only happens due to length
        assert_equal(
            len(sub[-1].trials[0].trajectory),
            len(traj_bad_extension)+self.dyn.n_frames_max-1
        )


        self.dyn.options['on_max_length'] = 'fail'


class TestAbstract(object):
    @raises_with_message_like(TypeError, "Can't instantiate abstract class")
    def test_abstract_pathmover(self):
        mover = paths.PathMover()

    @raises_with_message_like(TypeError, "Can't instantiate abstract class")
    def test_abstract_samplemover(self):
        mover = paths.SampleMover()

    @raises_with_message_like(TypeError, "Can't instantiate abstract class")
    def test_abstract_enginemover(self):
        mover = paths.EngineMover()

    @raises_with_message_like(TypeError, "Can't instantiate abstract class")
    def test_abstract_selectionmover(self):
        mover = paths.SelectionMover()

    @raises_with_message_like(TypeError, "Can't instantiate abstract class")
    def test_abstract_subtrajectoryselectmover(self):
        mover = paths.SubtrajectorySelectMover()
