from openpathsampling.tests.test_helpers import (data_filename)
import openpathsampling as paths
import openpathsampling.engines.toy as toys
from openpathsampling.pathsimulators.sshooting_simulator import SShootingSimulation
import numpy as np
import os

class TestSShootingSimulation(object):
    def setup(self):
        # PES is one-dimensional zero function (y(x) = 0)
        pes = toys.LinearSlope(m=[0.0], c=[0.0])
        topology = toys.Topology(n_spatial=1, masses=[1.0], pes=pes)
        integrator = toys.LeapfrogVerletIntegrator(0.02)
        options = {
            'integ' : integrator,
            'n_frames_max' : 1000,
            'n_steps_per_frame' : 1
        }
        self.engine = toys.Engine(options=options, topology=topology)
        # test uses snapshots with different velocities
        self.initial_snapshots = [toys.Snapshot(
                                      coordinates=np.array([[0.0]]),
                                      velocities=np.array([[1.0]]),
                                      engine=self.engine),
                                  toys.Snapshot(
                                      coordinates=np.array([[0.1]]),
                                      velocities=np.array([[-1.0]]),
                                      engine=self.engine),
                                  toys.Snapshot(
                                      coordinates=np.array([[-0.2]]),
                                      velocities=np.array([[2.0]]),
                                      engine=self.engine)]
        # trajectory length is set to 100 steps
        self.l = 100
        # reaction coordinate is just x coordinate
        rc = paths.FunctionCV("Id", lambda snap : snap.coordinates[0][0])
        # state S: [-0.5, 0.5]
        self.state_S = paths.CVDefinedVolume(rc, -0.5, 0.5)
        # define state labels
        self.state_labels = {
            "S" : self.state_S,
            "NotS" :~self.state_S}
        # velocities are not randomized
        randomizer = paths.NoModification()

        self.filename = data_filename("sshooting_test.nc")
        self.storage = paths.Storage(self.filename, mode="w")
        self.storage.save(self.initial_snapshots)

        self.simulation = SShootingSimulation(
            storage=self.storage,
            engine=self.engine,
            state_S=self.state_S,
            randomizer=randomizer,
            initial_snapshots=self.initial_snapshots,
            trajectory_length=self.l
        )
        self.simulation.output_stream = open(os.devnull, 'w')

    def teardown(self):
        if os.path.isfile(self.filename):
            os.remove(self.filename)
        paths.EngineMover.default_engine = None

    def test_initialization(self):
        sim = self.simulation
        assert len(sim.initial_snapshots) == 3
        assert isinstance(sim.mover, paths.SequentialMover)

    def test_simulation_run(self):
        self.simulation.run(n_per_snapshot=1)
        assert len(self.simulation.storage.steps) == 3
        sim = self.simulation

        steps = sim.storage.steps
        # Check if trajectory length has not been altered.
        assert sim.trajectory_length == self.l
        # Check if state S has not been altered.
        assert sim.state_S == self.state_S
        for i in range(len(steps)):
            step = steps[i]
            # there should be two trials
            assert len(step.change.trials) == 2
            # backward shot (first trial) should have length l+1
            assert len(step.change.trials[0].trajectory) == self.l+1
            # total trajectory (second trial = backward + forward shot should
            # have length 2*l+1
            assert len(step.change.trials[1].trajectory) == 2*self.l+1
            # middle point of trajectory should be initial snapshot
            t = step.change.trials[1].trajectory
            np.testing.assert_almost_equal(
                t[self.l].coordinates, self.initial_snapshots[i].coordinates
            )
            # all trajectories have intermediate points in S
            traj_summary = t.summarize_by_volumes(self.state_labels)
            assert traj_summary[0][0] == 'NotS'
            assert traj_summary[1][0] == 'S'
            assert traj_summary[2][0] == 'NotS'
            # last mover should be forward mover of simulation
            expected_mover = self.simulation.forward_mover
            assert step.change.canonical.mover == expected_mover
