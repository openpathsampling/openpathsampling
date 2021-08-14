from __future__ import division
from __future__ import absolute_import
from builtins import str
from builtins import range
from past.utils import old_div
from builtins import object
from .test_helpers import (raises_with_message_like, data_filename,
                           CalvinistDynamics, make_1d_traj,
                           assert_items_equal)
from nose.tools import (assert_equal, assert_not_equal, raises,
                        assert_almost_equal, assert_true, assert_greater)
# from nose.plugins.skip import SkipTest

from openpathsampling.pathsimulators import *
import openpathsampling as paths
import openpathsampling.engines.toy as toys
import numpy as np
import os

import logging
logging.getLogger('openpathsampling.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.netcdfplus').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.ensemble').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.engines').setLevel(logging.CRITICAL)


class TestAbstract(object):
    @raises_with_message_like(TypeError, "Can't instantiate abstract class")
    def test_abstract_volume(self):
        PathSimulator()


class TestFullBootstrapping(object):
    def setup(self):
        paths.InterfaceSet._reset()
        self.cv = paths.FunctionCV("Id", lambda snap: snap.xyz[0][0])
        cv_neg = paths.FunctionCV("Neg", lambda snap: -snap.xyz[0][0])
        self.stateA = paths.CVDefinedVolume(self.cv, -1.0, 0.0)
        self.stateB = paths.CVDefinedVolume(self.cv, 1.0, 2.0)
        self.stateC = paths.CVDefinedVolume(self.cv, 3.0, 4.0)
        interfacesAB = paths.VolumeInterfaceSet(
            self.cv, -1.0, [0.0, 0.2, 0.4]
        )
        interfacesBC = paths.VolumeInterfaceSet(
            self.cv, 1.0, [2.0, 2.2, 2.4]
        )
        interfacesBA = paths.VolumeInterfaceSet(
            cv_neg, -1.0, [-1.0, -0.8, -0.6]
        )

        network = paths.MISTISNetwork([
            (self.stateA, interfacesAB, self.stateB),
            (self.stateB, interfacesBC, self.stateC),
            (self.stateB, interfacesBA, self.stateA)
        ])
        self.tisAB = network.input_transitions[(self.stateA, self.stateB)]
        self.tisBC = network.input_transitions[(self.stateB, self.stateC)]
        self.tisBA = network.input_transitions[(self.stateB, self.stateA)]
        self.network = network
        self.snapA = make_1d_traj([-0.5])[0]

        self.noforbid_noextra_AB = paths.FullBootstrapping(
            transition=self.tisAB,
            snapshot=self.snapA
        )

    @raises(RuntimeError)
    def test_initial_max_length(self):
        engine = CalvinistDynamics([-0.5, -0.4, -0.3, -0.2, -0.1, 0.1, -0.1])
        bootstrap_AB_maxlength = paths.FullBootstrapping(
            transition=self.tisAB,
            snapshot=self.snapA,
            initial_max_length=3,
            engine=engine
        )
        bootstrap_AB_maxlength.output_stream = open(os.devnull, "w")
        bootstrap_AB_maxlength.run(build_attempts=1)

    def test_first_traj_ensemble(self):
        traj_starts_in = make_1d_traj([-0.2, -0.1, 0.1, -0.1])
        traj_starts_out = make_1d_traj([0.1, -0.1, 0.1, -0.1])
        traj_not_good = make_1d_traj([0.1, -0.1, 0.1])
        first_traj_ens = self.noforbid_noextra_AB.first_traj_ensemble
        assert_equal(first_traj_ens(traj_starts_in), True)
        assert_equal(first_traj_ens(traj_starts_out), True)
        assert_equal(first_traj_ens(traj_not_good), False)

    def test_sampling_ensembles(self):
        traj1 = make_1d_traj([-0.2, -0.1, 0.1, -0.1])
        traj2 = make_1d_traj([-0.1, 0.1, -0.1])
        traj3 = make_1d_traj([-0.1, 0.1, 0.3, -0.1])
        traj4 = make_1d_traj([0.1, 0.3, 0.1])
        all_ensembles = self.noforbid_noextra_AB.all_ensembles
        assert_equal(len(all_ensembles), 3)
        for ens in all_ensembles:
            assert_equal(ens(traj1), False)
            assert_equal(ens(traj4), False)
        assert_equal(all_ensembles[0](traj2), True)
        assert_equal(all_ensembles[0](traj3), True)
        assert_equal(all_ensembles[1](traj2), False)
        assert_equal(all_ensembles[1](traj3), True)
        assert_equal(all_ensembles[2](traj2), False)
        assert_equal(all_ensembles[2](traj3), False)

    def test_run_already_satisfied(self):
        engine = CalvinistDynamics([-0.5, 0.8, -0.1])
        bootstrap = FullBootstrapping(
            transition=self.tisAB,
            snapshot=self.snapA,
            engine=engine
        )
        bootstrap.output_stream = open(os.devnull, "w")
        gs = bootstrap.run()
        assert_equal(len(gs), 3)

    def test_run_extra_interfaces(self):
        engine = CalvinistDynamics([-0.5, 0.8, -0.1])
        bootstrap = FullBootstrapping(
            transition=self.tisAB,
            snapshot=self.snapA,
            engine=engine,
            extra_interfaces=[paths.CVDefinedVolume(self.cv, -1.0, 0.6)]
        )
        bootstrap.output_stream = open(os.devnull, "w")
        gs = bootstrap.run()
        assert_equal(len(gs), 4)

    def test_run_forbidden_states(self):
        engine = CalvinistDynamics([-0.5, 0.3, 3.2, -0.1, 0.8, -0.1])
        # first, without setting forbidden_states
        bootstrap1 = FullBootstrapping(
            transition=self.tisAB,
            snapshot=self.snapA,
            engine=engine
        )
        bootstrap1.output_stream = open(os.devnull, "w")
        gs1 = bootstrap1.run()
        assert_equal(len(gs1), 3)
        assert_items_equal(self.cv(gs1[0]), [-0.5, 0.3, 3.2, -0.1])
        # now with setting forbidden_states
        bootstrap2 = FullBootstrapping(
            transition=self.tisAB,
            snapshot=self.snapA,
            engine=engine,
            forbidden_states=[self.stateC]
        )
        bootstrap2.output_stream = open(os.devnull, "w")
        # make sure this is where we get the error
        try:
            bootstrap2.run()
        except RuntimeError:
            pass

    @raises(RuntimeError)
    def test_too_much_bootstrapping(self):
        engine = CalvinistDynamics([-0.5, 0.2, -0.1])
        bootstrap = FullBootstrapping(
            transition=self.tisAB,
            snapshot=self.snapA,
            engine=engine,
        )
        bootstrap.output_stream = open(os.devnull, "w")
        bootstrap.run(max_ensemble_rounds=1)


class TestShootFromSnapshotsSimulation(object):
    # note that most of ShootFromSnapshotSimulation is tested in the tests
    # for CommittorSimulation. This is just an additional test to show that
    # using different ensembles from the ones used for the committor will
    # also work.
    def setup(self):
        # As a test system, let's use 1D motion on a flat potential. If the
        # velocity is positive, you right the state on the right. If it is
        # negative, you hit the state on the left.
        pes = toys.LinearSlope(m=[0.0], c=[0.0])  # flat line
        topology = toys.Topology(n_spatial=1, masses=[1.0], pes=pes)
        integrator = toys.LeapfrogVerletIntegrator(0.1)
        options = {
            'integ': integrator,
            'n_frames_max': 100000,
            'n_steps_per_frame': 5
        }
        self.engine = toys.Engine(options=options, topology=topology)
        self.snap0 = toys.Snapshot(coordinates=np.array([[0.0]]),
                                   velocities=np.array([[1.0]]),
                                   engine=self.engine)
        cv = paths.FunctionCV("Id", lambda snap: snap.coordinates[0][0])
        starting_volume = paths.CVDefinedVolume(cv, -0.01, 0.01)
        forward_ensemble = paths.LengthEnsemble(5)
        backward_ensemble = paths.LengthEnsemble(3)
        randomizer = paths.NoModification()

        self.filename = data_filename("shoot_from_snaps.nc")
        self.storage = paths.Storage(self.filename, 'w')
        self.simulation = ShootFromSnapshotsSimulation(
            storage=self.storage,
            engine=self.engine,
            starting_volume=starting_volume,
            forward_ensemble=forward_ensemble,
            backward_ensemble=backward_ensemble,
            randomizer=randomizer,
            initial_snapshots=self.snap0
        )
        self.simulation.output_stream = open(os.devnull, "w")

    def teardown(self):
        if os.path.isfile(self.filename):
            os.remove(self.filename)
        paths.EngineMover.default_engine = None

    def test_run_arbitrary_ensemble(self):
        # integration test of the whole thing, including storage
        self.simulation.run(10)
        self.storage.close()
        analysis = paths.Storage(self.filename, 'r')
        sim = analysis.pathsimulators[0]
        assert_equal(len(analysis.steps), 10)
        length_to_submover = {5: [], 3: []}
        for step in analysis.steps:
            step.active.sanity_check()
            assert_equal(len(step.active), 1)
            active_sample = step.active[0]
            change = step.change
            assert_equal(change.mover, sim.mover)
            # KeyError here indicates problem with lengths generated
            length_to_submover[len(active_sample)] += change.subchange.mover

        for k in length_to_submover:
            # allow 0 or 1  because maybe we made no trials with submover
            assert_true(len(set(length_to_submover[k])) <= 1)


class TestCommittorSimulation(object):
    def setup(self):
        # As a test system, let's use 1D motion on a flat potential. If the
        # velocity is positive, you right the state on the right. If it is
        # negative, you hit the state on the left.
        pes = toys.LinearSlope(m=[0.0], c=[0.0])  # flat line
        topology = toys.Topology(n_spatial=1, masses=[1.0], pes=pes)
        integrator = toys.LeapfrogVerletIntegrator(0.1)
        options = {
            'integ': integrator,
            'n_frames_max': 100000,
            'n_steps_per_frame': 5
        }
        self.engine = toys.Engine(options=options, topology=topology)
        self.snap0 = toys.Snapshot(coordinates=np.array([[0.0]]),
                                   velocities=np.array([[1.0]]),
                                   engine=self.engine)
        cv = paths.FunctionCV("Id", lambda snap: snap.coordinates[0][0])
        self.left = paths.CVDefinedVolume(cv, float("-inf"), -1.0)
        self.right = paths.CVDefinedVolume(cv, 1.0, float("inf"))
        self.state_labels = {"Left": self.left,
                             "Right": self.right,
                             "None": ~(self.left | self.right)}

        randomizer = paths.NoModification()

        self.filename = data_filename("committor_test.nc")
        self.storage = paths.Storage(self.filename, mode="w")
        self.storage.save(self.snap0)

        self.simulation = CommittorSimulation(storage=self.storage,
                                              engine=self.engine,
                                              states=[self.left, self.right],
                                              randomizer=randomizer,
                                              initial_snapshots=self.snap0)
        self.simulation.output_stream = open(os.devnull, 'w')

    def teardown(self):
        self.storage.close()
        if os.path.isfile(self.filename):
            os.remove(self.filename)
        paths.EngineMover.default_engine = None

    def test_initialization(self):
        sim = self.simulation  # convenience
        assert_equal(len(sim.initial_snapshots), 1)
        assert_true(isinstance(sim.mover, paths.RandomChoiceMover))

    def test_storage(self):
        self.storage.tag['simulation'] = self.simulation
        self.storage.close()
        read_store = paths.Storage(self.filename, 'r')
        sim = read_store.tag['simulation']
        new_filename = data_filename("test2.nc")
        sim.storage = paths.Storage(new_filename, 'w')
        sim.output_stream = open(os.devnull, 'w')
        sim.run(n_per_snapshot=2)
        sim.storage.close()
        if os.path.isfile(new_filename):
            os.remove(new_filename)
        self.storage = read_store  # teardown will get rid of this

    def test_committor_run(self):
        self.simulation.run(n_per_snapshot=20)
        assert_equal(len(self.simulation.storage.steps), 20)
        counts = {'fwd': 0, 'bkwd': 0}
        for step in self.simulation.storage.steps:
            step.active.sanity_check()  # traj is in ensemble
            traj = step.active[0].trajectory
            traj_str = traj.summarize_by_volumes_str(self.state_labels)
            if traj_str == "None-Right":
                assert_equal(step.change.canonical.mover,
                             self.simulation.forward_mover)
                assert_equal(step.active[0].ensemble,
                             self.simulation.forward_ensemble)
                counts['fwd'] += 1
            elif traj_str == "Left-None":
                assert_equal(step.change.canonical.mover,
                             self.simulation.backward_mover)
                assert_equal(step.active[0].ensemble,
                             self.simulation.backward_ensemble)
                counts['bkwd'] += 1
            else:
                raise AssertionError(
                    str(traj_str) + "is neither 'None-Right' nor 'Left-None'"
                )
        assert_true(counts['fwd'] > 0)
        assert_true(counts['bkwd'] > 0)
        assert_equal(counts['fwd'] + counts['bkwd'], 20)

    def test_forward_only_committor(self):
        sim = CommittorSimulation(storage=self.storage,
                                  engine=self.engine,
                                  states=[self.left, self.right],
                                  randomizer=paths.NoModification(),
                                  initial_snapshots=self.snap0,
                                  direction=1)
        sim.output_stream = open(os.devnull, 'w')
        sim.run(n_per_snapshot=10)
        assert_equal(len(sim.storage.steps), 10)
        for step in self.simulation.storage.steps:
            s = step.active[0]
            step.active.sanity_check()  # traj is in ensemble
            assert_equal(
                s.trajectory.summarize_by_volumes_str(self.state_labels),
                "None-Right"
            )
            assert_equal(s.ensemble, sim.forward_ensemble)
            assert_equal(step.change.canonical.mover,
                         sim.forward_mover)

    def test_backward_only_committor(self):
        sim = CommittorSimulation(storage=self.storage,
                                  engine=self.engine,
                                  states=[self.left, self.right],
                                  randomizer=paths.NoModification(),
                                  initial_snapshots=self.snap0,
                                  direction=-1)
        sim.output_stream = open(os.devnull, 'w')
        sim.run(n_per_snapshot=10)
        assert_equal(len(sim.storage.steps), 10)
        for step in self.simulation.storage.steps:
            s = step.active[0]
            step.active.sanity_check()  # traj is in ensemble
            assert_equal(
                s.trajectory.summarize_by_volumes_str(self.state_labels),
                "Left-None"
            )
            assert_equal(s.ensemble, sim.backward_ensemble)
            assert_equal(step.change.canonical.mover,
                         sim.backward_mover)

    def test_multiple_initial_snapshots(self):
        snap1 = toys.Snapshot(coordinates=np.array([[0.1]]),
                              velocities=np.array([[-1.0]]),
                              engine=self.engine)
        sim = CommittorSimulation(storage=self.storage,
                                  engine=self.engine,
                                  states=[self.left, self.right],
                                  randomizer=paths.NoModification(),
                                  initial_snapshots=[self.snap0, snap1])
        sim.output_stream = open(os.devnull, 'w')
        sim.run(10)
        assert_equal(len(self.storage.steps), 20)
        snap0_coords = self.snap0.coordinates.tolist()
        snap1_coords = snap1.coordinates.tolist()
        count = {self.snap0: 0, snap1: 0}
        for step in self.storage.steps:
            shooting_snap = step.change.canonical.details.shooting_snapshot
            if shooting_snap.coordinates.tolist() == snap0_coords:
                mysnap = self.snap0
            elif shooting_snap.coordinates.tolist() == snap1_coords:
                mysnap = snap1
            else:
                msg = "Shooting snapshot matches neither test snapshot"
                raise AssertionError(msg)
            count[mysnap] += 1
        assert_equal(count, {self.snap0: 10, snap1: 10})

    def test_randomized_committor(self):
        # this shows that we get both states even with forward-only
        # shooting, if the randomizer gives the negative velocities
        randomizer = paths.RandomVelocities(beta=1.0)
        sim = CommittorSimulation(storage=self.storage,
                                  engine=self.engine,
                                  states=[self.left, self.right],
                                  randomizer=randomizer,
                                  initial_snapshots=self.snap0,
                                  direction=1)
        sim.output_stream = open(os.devnull, 'w')
        sim.run(50)
        assert_equal(len(sim.storage.steps), 50)
        counts = {'None-Right': 0,
                  'Left-None': 0,
                  'None-Left': 0,
                  'Right-None': 0}
        for step in sim.storage.steps:
            step.active.sanity_check()  # traj is in ensemble
            traj = step.active[0].trajectory
            traj_str = traj.summarize_by_volumes_str(self.state_labels)
            try:
                counts[traj_str] += 1
            except KeyError:
                msg = "Got trajectory described as '{0}', length {1}"
                # this might be okay if it is 'None', length 100000
                raise AssertionError(msg.format(traj_str, len(traj)))
        assert_equal(counts['Left-None'], 0)
        assert_equal(counts['Right-None'], 0)
        assert_true(counts['None-Left'] > 0)
        assert_true(counts['None-Right'] > 0)
        assert_equal(sum(counts.values()), 50)


class TestReactiveFluxSimulation(object):
    def setup(self):
        # PES is one-dimensional linear slope (y(x) = x)
        pes = toys.LinearSlope(m=[-1.0], c=[0.0])
        # one particle with mass 1.0
        topology = toys.Topology(n_spatial=1, masses=[1.0], pes=pes)
        integrator = toys.LeapfrogVerletIntegrator(0.02)
        options = {
            'integ' : integrator,
            'n_frames_max' : 1000,
            'n_steps_per_frame' : 5
        }
        self.engine = toys.Engine(options=options, topology=topology)
        # test uses three snapshots with different velocities
        # 0: direction ok, velocity too low => falls back to dividing surface
        # 1: wrong direction => backward shot towards B
        # 2: direction ok, velocity high enough => successfull new trajectory
        self.initial_snapshots = [toys.Snapshot(
                                      coordinates=np.array([[0.0]]),
                                      velocities=np.array([[1.0]]),
                                      engine=self.engine),
                                  toys.Snapshot(
                                      coordinates=np.array([[0.0]]),
                                      velocities=np.array([[-1.0]]),
                                      engine=self.engine),
                                  toys.Snapshot(
                                      coordinates=np.array([[0.0]]),
                                      velocities=np.array([[2.0]]),
                                      engine=self.engine)]
        # reaction coordinate is just x coordinate
        rc = paths.FunctionCV("Id", lambda snap : snap.coordinates[0][0])
        # state A: [-inf, -1]
        self.state_A = paths.CVDefinedVolume(rc, float("-inf"), -1.0)
        # area between A and dividing surface: [-1, 0]
        self.towards_A = paths.CVDefinedVolume(rc, -1.0, 0.0)
        # state B: [1, inf]
        self.state_B = paths.CVDefinedVolume(rc, 1.0, float("inf"))
        # define state labels
        self.state_labels = {
            "A" : self.state_A,
            "B" : self.state_B,
            "ToA": self.towards_A,
            "None" :~(self.state_A | self.state_B | self.towards_A)}

        # velocities are not randomized
        randomizer = paths.NoModification()

        self.filename = data_filename("rf_test.nc")
        self.storage = paths.Storage(self.filename, mode="w")
        self.storage.save(self.initial_snapshots)

        self.simulation = ReactiveFluxSimulation(
                              storage=self.storage,
                              engine=self.engine,
                              states=[self.state_A, self.state_B],
                              randomizer=randomizer,
                              initial_snapshots=self.initial_snapshots,
                              rc=rc)
        self.simulation.output_stream = open(os.devnull, 'w')

    def teardown(self):
        if os.path.isfile(self.filename):
            os.remove(self.filename)
        paths.EngineMover.default_engine = None

    def test_initialization(self):
        sim = self.simulation
        assert_equal(len(sim.initial_snapshots), 3)
        assert_true(isinstance(sim.mover, paths.ConditionalSequentialMover))

    def test_simulation_run(self):
        self.simulation.run(n_per_snapshot=1)
        assert_equal(len(self.simulation.storage.steps), 3)

        # snapshot 0, fails at backward shot (falls back to dividing surface)
        step = self.simulation.storage.steps[0]
        # last mover should be backward_mover of simulation
        assert_equal(step.change.canonical.mover,
                     self.simulation.backward_mover)
        # active ensemble should be starting ensemble
        assert_equal(step.active[0].ensemble,
                     self.simulation.starting_ensemble)
        # analyze trajectory, last step should be in 'None', the rest in 'ToA'
        traj = step.change.trials[0].trajectory
        traj_summary = traj.summarize_by_volumes(self.state_labels)
        assert_equal(traj_summary[0], ('None', 1))
        assert_equal(traj_summary[1][0], 'ToA')
        assert_greater(traj_summary[1][1], 1)

        # snapshot 1, fails at backward shot (wrong direction)
        step = self.simulation.storage.steps[1]
        # last mover should be backward_mover of simulation
        assert_equal(step.change.canonical.mover,
                     self.simulation.backward_mover)
        # active ensemble should be starting ensemble
        assert_equal(step.active[0].ensemble,
                     self.simulation.starting_ensemble)
        # analyze trajectory, backwards trajectory reaches immediately 'None'
        traj = step.change.trials[0].trajectory
        traj_summary = traj.summarize_by_volumes(self.state_labels)
        assert_equal(traj_summary[0], ('None', 2))

        # snapshot 2, is accepted
        step = self.simulation.storage.steps[2]
        # last mover should be forward_mover of simulation
        assert_equal(step.change.canonical.mover,
                     self.simulation.forward_mover)
        # active ensemble should not be starting ensemble
        assert_not_equal(step.active[0].ensemble,
                     self.simulation.starting_ensemble)
        # analyze active trajectory, trajectory should start in 'A', end in 'B'
        traj = step.active[0].trajectory
        traj_summary = traj.summarize_by_volumes(self.state_labels)
        assert_equal(traj_summary[0], ('A', 1))
        assert_equal(traj_summary[1][0], 'ToA')
        assert_greater(traj_summary[1][1], 1)
        assert_equal(traj_summary[2][0], 'None')
        assert_greater(traj_summary[2][1], 1)
        assert_equal(traj_summary[3], ('B', 1))


class TestDirectSimulation(object):
    def setup(self):
        pes = toys.HarmonicOscillator(A=[1.0], omega=[1.0], x0=[0.0])
        topology = toys.Topology(n_spatial=1, masses=[1.0], pes=pes)
        integrator = toys.LeapfrogVerletIntegrator(0.1)
        options = {
            'integ': integrator,
            'n_frames_max': 100000,
            'n_steps_per_frame': 2
        }
        self.engine = toys.Engine(options=options, topology=topology)
        self.snap0 = toys.Snapshot(coordinates=np.array([[0.0]]),
                                   velocities=np.array([[1.0]]),
                                   engine=self.engine)
        cv = paths.FunctionCV("Id", lambda snap: snap.coordinates[0][0])
        self.cv = cv
        self.center = paths.CVDefinedVolume(cv, -0.2, 0.2)
        self.interface = paths.CVDefinedVolume(cv, -0.3, 0.3)
        self.outside = paths.CVDefinedVolume(cv, 0.6, 0.9)
        self.extra = paths.CVDefinedVolume(cv, -1.5, -0.9)
        self.flux_pairs = [(self.center, self.interface)]
        self.sim = DirectSimulation(storage=None,
                                    engine=self.engine,
                                    states=[self.center, self.outside],
                                    flux_pairs=self.flux_pairs,
                                    initial_snapshot=self.snap0)

    def test_run(self):
        self.sim.run(200)
        assert_true(len(self.sim.transition_count) > 1)
        assert_true(len(self.sim.flux_events[self.flux_pairs[0]]) > 1)

    def test_results(self):
        self.sim.run(200)
        results = self.sim.results
        assert_equal(len(results), 2)
        assert_equal(set(results.keys()),
                     set(['transition_count', 'flux_events']))
        assert_equal(results['transition_count'], self.sim.transition_count)
        assert_equal(results['flux_events'], self.sim.flux_events)

    def test_load_results(self):
        left_interface = paths.CVDefinedVolume(self.cv, -0.3, float("inf"))
        right_interface = paths.CVDefinedVolume(self.cv, float("-inf"), 0.3)
        fake_transition_count = [
            (self.center, 1), (self.outside, 4), (self.center, 7),
            (self.extra, 10), (self.center, 12), (self.outside, 14)
        ]
        fake_flux_events = {(self.center, right_interface):
                            [(15, 3), (23, 15), (48, 23)],
                            (self.center, left_interface):
                            [(97, 34), (160, 97)]}
        results = {'transition_count': fake_transition_count,
                   'flux_events': fake_flux_events}
        self.sim.load_results(results)
        assert_equal(self.sim.transition_count, fake_transition_count)
        assert_equal(self.sim.flux_events, fake_flux_events)

    def test_transitions(self):
        # set fake data
        self.sim.transition_count = [
            (self.center, 1), (self.outside, 4), (self.center, 7),
            (self.extra, 10), (self.center, 12), (self.outside, 14)
        ]
        assert_equal(self.sim.n_transitions,
                     {(self.center, self.outside): 2,
                      (self.outside, self.center): 1,
                      (self.center, self.extra): 1,
                      (self.extra, self.center): 1})
        assert_equal(self.sim.transitions,
                     {(self.center, self.outside): [3, 2],
                      (self.outside, self.center): [3],
                      (self.center, self.extra): [3],
                      (self.extra, self.center): [2]})

    def test_rate_matrix(self):
        self.sim.states += [self.extra]
        self.sim.transition_count = [
            (self.center, 1), (self.outside, 4), (self.center, 7),
            (self.extra, 10), (self.center, 12), (self.outside, 14)
        ]
        step_size = self.sim.engine.snapshot_timestep
        # As of pandas 0.18.1, callables can be used in `df.loc`, etc. Since
        # we're using (callable) volumes for labels of columns/indices in
        # our dataframes, this sucks for us. Raise an issue with pandas?
        rate_matrix = self.sim.rate_matrix.values
        # order is center, outside, extra
        nan = float("nan")
        t_center = 3.0 + 3.0 + 2.0
        t_outside = 3.0
        t_extra = 2.0
        test_matrix = np.array([[nan, 2.0/t_center, 1.0/t_center],
                                [1.0/t_outside, nan, nan],
                                [1.0/t_extra, nan, nan]]) / step_size
        # for some reason, np.testing.assert_allclose(..., equal_nan=True)
        # was raising errors on this input. this hack gets the behavior
        for i in range(len(self.sim.states)):
            for j in range(len(self.sim.states)):
                if np.isnan(test_matrix[i][j]):
                    assert_true(np.isnan(rate_matrix[i][j]))
                else:
                    assert_almost_equal(rate_matrix[i][j],
                                        test_matrix[i][j])

    def test_fluxes(self):
        left_interface = paths.CVDefinedVolume(self.cv, -0.3, float("inf"))
        right_interface = paths.CVDefinedVolume(self.cv, float("-inf"), 0.3)
        sim = DirectSimulation(storage=None,
                               engine=self.engine,
                               states=[self.center, self.outside],
                               flux_pairs=[(self.center, left_interface),
                                           (self.center, right_interface)],
                               initial_snapshot=self.snap0)
        fake_flux_events = {(self.center, right_interface):
                            [(15, 3), (23, 15), (48, 23)],
                            (self.center, left_interface):
                            [(97, 34), (160, 97)]}
        sim.flux_events = fake_flux_events
        n_flux_events = {(self.center, right_interface): 3,
                         (self.center, left_interface): 2}
        assert_equal(sim.n_flux_events, n_flux_events)
        time_step = sim.engine.snapshot_timestep
        expected_fluxes = {(self.center, right_interface):
                           1.0 / (((15-3) + (23-15) + (48-23))/3.0) / time_step,
                           (self.center, left_interface):
                           1.0 / (((97-34) + (160-97))/2.0) / time_step}
        for p in expected_fluxes:
            assert_almost_equal(sim.fluxes[p], expected_fluxes[p])

    def test_multiple_cv_flux(self):
        # To check for the multiple interface set case, we need to have two
        # dimensions. We can hack two "independent" dimensions from a one
        # dimensional system by making the second CV non-monotonic with the
        # first. For the full trajectory, we need snapshots `S` (in the
        # state); `I` (interstitial: outside the state, but not outside
        # either interface); `X_a` (outside interface alpha, not outside
        # interface beta); `X_b` (outside interface beta, not outside
        # interface alpha); and `X_ab` (outside interface alpha and beta).
        cv1 = self.cv
        cv2 = paths.FunctionCV("abs_sin",
                               lambda snap: np.abs(np.sin(snap.xyz[0][0])))
        state = paths.CVDefinedVolume(cv1,
                                      old_div(-np.pi, 8.0),
                                      old_div(np.pi, 8.0))
        other_state = paths.CVDefinedVolume(cv1,
                                            -5.0/8.0 * np.pi,
                                            -3.0/8.0 * np.pi)
        alpha = paths.CVDefinedVolume(cv1, float("-inf"), 3.0/8.0*np.pi)
        beta = paths.CVDefinedVolume(cv2,
                                     float("-inf"),
                                     old_div(np.sqrt(2), 2.0))
        # approx     alpha: x < 1.17   beta: abs(sin(x)) < 0.70
        S = 0              # cv1 =  0.00; cv2 = 0.00
        I = old_div(np.pi, 5.0)      # cv1 =  0.63; cv2 = 0.59
        X_a = np.pi        # cv1 =  3.14; cv2 = 0.00
        X_b = old_div(-np.pi, 3.0)    # cv1 = -1.05; cv2 = 0.87
        X_ab = old_div(np.pi, 2.0)    # cv1 =  1.57; cv2 = 1.00
        other = old_div(-np.pi, 2.0)  # cv1 = -1.57; cv2 = 1.00
        # That hack is utterly crazy, but I'm kinda proud of it!
        predetermined = [S, S, I, X_a,   # (2) first exit
                         S, X_a,         # (4) cross A
                         S, X_ab,        # (6) cross A & B
                         I, S, X_b,      # (9) cross B
                         S, I, X_b,      # (12) cross B
                         other, I, X_b,  # (15) cross to other state
                         S, X_b,         # (17) first cross B
                         S, X_a,         # (19) first cross A
                         S, S, X_ab,     # (22) cross A & B
                         I, X_ab,        # (24) recrossing test
                         S, I,           # (26) false crossing test
                         S, S]
        engine = CalvinistDynamics(predetermined)
        init = make_1d_traj([S])
        sim = DirectSimulation(storage=None,
                               engine=engine,
                               states=[state, other_state],
                               flux_pairs=[(state, alpha), (state, beta)],
                               initial_snapshot=init[0])
        sim.run(len(predetermined)-1)
        # subtract 1 from the indices in `predetermined`, b/c 0 index of the
        # traj comes after the found initial step
        expected_flux_events = {
            (state, alpha): [(4, 2), (6, 4), (22, 19)],
            (state, beta): [(9, 6), (12, 9), (22, 17)]
        }
        assert_equal(len(sim.flux_events), 2)
        assert_equal(sim.flux_events[(state, alpha)],
                     expected_flux_events[(state, alpha)])
        assert_equal(sim.flux_events[(state, beta)],
                     expected_flux_events[(state, beta)])

    def test_simple_flux(self):
        state = self.center
        interface = self.interface
        A = 0.0
        I = 0.25
        X = 0.35
        # 0 index of traj comes after first found step in CalvinistDyn
        # FRAME INDEX:      0              5             10    12
        predetermined = [X, A, I, X, I, X, I, A, I, A, I, X, I, A]
        engine = CalvinistDynamics(predetermined)
        init = make_1d_traj([X])
        sim = DirectSimulation(storage=None,
                               engine=engine,
                               states=[state],
                               flux_pairs=[(state, interface)],
                               initial_snapshot=init[0])
        sim.run(len(predetermined) - 1)
        expected_flux_events = {(state, interface): [(10, 2)]}
        assert_equal(len(sim.flux_events), 1)
        assert_equal(sim.flux_events[(state, interface)],
                     expected_flux_events[(state, interface)])

    def test_sim_with_storage(self):
        tmpfile = data_filename("direct_sim_test.nc")
        if os.path.isfile(tmpfile):
            os.remove(tmpfile)

        storage = paths.Storage(tmpfile, "w", self.snap0)
        sim = DirectSimulation(storage=storage,
                               engine=self.engine,
                               states=[self.center, self.outside],
                               initial_snapshot=self.snap0)

        sim.run(200)
        storage.close()
        read_store = paths.AnalysisStorage(tmpfile)
        assert_equal(len(read_store.trajectories), 1)
        traj = read_store.trajectories[0]
        assert_equal(len(traj), 201)
        read_store.close()
        os.remove(tmpfile)


class TestPathSampling(object):
    def setup(self):
        paths.InterfaceSet._reset()
        self.cv = paths.FunctionCV("x", lambda x: x.xyz[0][0])
        self.state_A = paths.CVDefinedVolume(self.cv, float("-inf"), 0.0)
        self.state_B = paths.CVDefinedVolume(self.cv, 1.0, float("inf"))
        pes = paths.engines.toy.LinearSlope([0, 0, 0], 0)
        integ = paths.engines.toy.LangevinBAOABIntegrator(0.01, 0.1, 2.5)
        topology = paths.engines.toy.Topology(n_spatial=3, masses=[1.0],
                                              pes=pes)
        self.engine = paths.engines.toy.Engine(options={'integ': integ},
                                               topology=topology)

        interfaces = paths.VolumeInterfaceSet(self.cv, float("-inf"),
                                              [0.0, 0.1, 0.2])
        network = paths.MISTISNetwork([
            (self.state_A, interfaces, self.state_B)
        ])
        init_traj = make_1d_traj([-0.1, 0.2, 0.5, 0.8, 1.1])
        scheme = paths.MoveScheme(network)
        scheme.append([
            paths.strategies.OneWayShootingStrategy(
                selector=paths.UniformSelector(),
                engine=self.engine
            ),
            paths.strategies.PathReversalStrategy(),
            paths.strategies.OrganizeByMoveGroupStrategy()
        ])
        init_cond = scheme.initial_conditions_from_trajectories(init_traj)
        self.scheme = scheme
        self.init_cond = init_cond
        self.sim = PathSampling(storage=None, move_scheme=scheme,
                                sample_set=init_cond)

    def test_run_until_decorrelated(self):
        def all_snaps(sample_set):
            return set(sum([s.trajectory for s in sample_set], []))
        initial_snaps = all_snaps(self.sim.sample_set)
        self.sim.run_until_decorrelated()
        final_snaps = all_snaps(self.sim.sample_set)
        assert initial_snaps & final_snaps == set([])
        # test time reversal
        init_xyz = set(s.xyz.tobytes() for s in initial_snaps)
        final_xyz = set(s.xyz.tobytes() for s in final_snaps)
        assert init_xyz & final_xyz == set([])

    def test_save_initial_scheme(self, tmpdir):
        # check that we actually save scheme when we save this
        filename = tmpdir.join("temp.nc")
        storage = paths.Storage(str(filename), mode='w')
        assert len(storage.schemes) == 0
        sim = paths.PathSampling(storage=storage,
                                 move_scheme=self.scheme,
                                 sample_set=self.init_cond)
        assert len(storage.schemes) == 1
