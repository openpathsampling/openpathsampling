from test_helpers import (raises_with_message_like, data_filename,
                          CalvinistDynamics, make_1d_traj)
from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        raises, assert_almost_equal, assert_true)
from nose.plugins.skip import SkipTest

from openpathsampling.pathsimulator import *
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

class testAbstract(object):
    @raises_with_message_like(TypeError, "Can't instantiate abstract class")
    def test_abstract_volume(self):
        mover = PathSimulator()


class testFullBootstrapping(object):
    def setup(self):
        self.cv = paths.CV_Function("Id", lambda snap: snap.xyz[0][0])
        cv_neg = paths.CV_Function("Neg", lambda snap: -snap.xyz[0][0])
        self.stateA = paths.CVRangeVolume(self.cv, -1.0, 0.0)
        self.stateB = paths.CVRangeVolume(self.cv, 1.0, 2.0)
        self.stateC = paths.CVRangeVolume(self.cv, 3.0, 4.0)
        interfacesAB = paths.VolumeFactory.CVRangeVolumeSet(
            self.cv, -1.0, [0.0, 0.2, 0.4]
        )
        interfacesBC = paths.VolumeFactory.CVRangeVolumeSet(
            self.cv, 1.0, [2.0, 2.2, 2.4]
        )
        interfacesBA = paths.VolumeFactory.CVRangeVolumeSet(
            cv_neg, -1.0, [-1.0, -0.8, -0.6]
        )

        network = paths.MISTISNetwork([
            (self.stateA, interfacesAB, self.cv, self.stateB),
            (self.stateB, interfacesBC, self.cv, self.stateC),
            (self.stateB, interfacesBA, cv_neg, self.stateA)
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
            initial_max_length = 3,
            engine=engine
        )
        bootstrap_AB_maxlength.output_stream = open(os.devnull, "w")
        gs = bootstrap_AB_maxlength.run(build_attempts=1)

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
        raise SkipTest

    def test_run_extra_interfaces(self):
        raise SkipTest

    def test_run_forbidden_states(self):
        raise SkipTest

    def test_too_much_bootstrapping(self):
        raise SkipTest

class testCommittorSimulation(object):
    def setup(self):
        # As a test system, let's use 1D motion on a flat potential. If the
        # velocity is positive, you right the state on the right. If it is
        # negative, you hit the state on the left.
        pes = toys.LinearSlope(m=[0.0], c=[0.0]) # flat line
        topology = toys.Topology(n_spatial=1, masses=[1.0], pes=pes)
        self.snap0 = toys.Snapshot(coordinates=np.array([[0.0]]),
                                   velocities=np.array([[1.0]]),
                                   topology=topology)
        integrator = toys.LeapfrogVerletIntegrator(0.1)
        options = {
            'integ': integrator,
            'n_frames_max': 100000,
            'nsteps_per_frame': 5
        }
        self.engine = toys.Engine(options=options, template=self.snap0)
        cv = paths.CV_Function("Id", lambda snap : snap.coordinates[0][0])
        self.left = paths.CVRangeVolume(cv, float("-inf"), -1.0)
        self.right = paths.CVRangeVolume(cv, 1.0, float("inf"))
        self.state_labels = {"Left" : self.left,
                             "Right" : self.right,
                             "None" : ~(self.left | self.right)}

        randomizer = paths.NoModification()

        self.filename = data_filename("committor_test.nc")
        self.storage = paths.Storage(self.filename, 
                                     mode="w", 
                                     template=self.snap0)

        self.simulation = CommittorSimulation(storage=self.storage,
                                              engine=self.engine,
                                              states=[self.left, self.right],
                                              randomizer=randomizer,
                                              initial_snapshots=self.snap0)

    def teardown(self):
        import os
        if os.path.isfile(self.filename):
            os.remove(self.filename)
        paths.EngineMover.default_engine = None

    def test_initialization(self):
        sim = self.simulation  # convenience
        assert_equal(len(sim.initial_snapshots), 1)
        assert_true(isinstance(sim.mover, paths.RandomChoiceMover))

    def test_committor_run(self):
        self.simulation.run(n_per_snapshot=20)
        assert_equal(len(self.simulation.storage.steps), 20)
        counts = {'fwd' : 0, 'bkwd' : 0}
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
                              topology=self.snap0.topology)
        sim = CommittorSimulation(storage=self.storage,
                                  engine=self.engine,
                                  states=[self.left, self.right],
                                  randomizer=paths.NoModification(),
                                  initial_snapshots=[self.snap0, snap1])
        sim.run(10)
        assert_equal(len(self.storage.steps), 20)
        snap0_coords = self.snap0.coordinates.tolist()
        snap1_coords = snap1.coordinates.tolist()
        count = {self.snap0: 0, snap1: 0}
        for step in self.storage.steps:
            # TODO: this should in step.change.canonical.details
            shooting_snap = step.change.trials[0].details.shooting_snapshot
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
        sim.run(50)
        assert_equal(len(sim.storage.steps), 50)
        counts = {'None-Right' : 0,
                  'Left-None' : 0,
                  'None-Left' : 0,
                  'Right-None' : 0}
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
