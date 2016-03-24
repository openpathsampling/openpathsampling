from test_helpers import raises_with_message_like, data_filename
from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        raises, assert_almost_equal, assert_true)
from nose.plugins.skip import SkipTest

from openpathsampling.pathsimulator import *
import openpathsampling as paths
import openpathsampling.engines.toy as toys
import numpy as np

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


class testCommittorSimulation(object):
    def setup(self):
        # As a test system, let's use 1D motion on a flat potential. If the
        # velocity is positive, you right the state on the right. If it is
        # negative, you hit the state on the left.
        pes = toys.LinearSlope(m=[0.0], c=[0.0]) # flat line
        topology = toys.Topology(n_spatial=1, masses=[1.0], pes=pes)
        snap0 = toys.Snapshot(coordinates=np.array([[0.0]]),
                              velocities=np.array([[1.0]]),
                              topology=topology)
        integrator = toys.LeapfrogVerletIntegrator(0.1)
        options = {
            'integ': integrator,
            'n_frames_max': 1000,
            'nsteps_per_frame': 5
        }
        engine = toys.Engine(options=options, template=snap0)
        cv = paths.CV_Function("Id", lambda snap : snap.coordinates[0][0])
        self.left = paths.CVRangeVolume(cv, float("-inf"), -1.0)
        self.right = paths.CVRangeVolume(cv, 1.0, float("inf"))

        randomizer = paths.NoModification()

        self.filename = data_filename("committor_test.nc")
        storage = paths.Storage(self.filename, mode="w", template=snap0)

        self.simulation = CommittorSimulation(storage=storage,
                                              engine=engine,
                                              states=[self.left, self.right],
                                              randomizer=randomizer,
                                              initial_snapshots=snap0)

    def teardown(self):
        import os
        if os.path.isfile(self.filename):
            os.remove(self.filename)

    def test_initialization(self):
        sim = self.simulation  # convenience
        assert_equal(len(sim.initial_snapshots), 1)
        assert_true(isinstance(sim.mover, paths.RandomChoiceMover))

    def test_committor_run(self):
        self.simulation.run(n_per_snapshot=10)
        assert_equal(len(self.simulation.storage.steps), 10)
        state_label = {"Left" : self.left,
                       "Right" : self.right, 
                       "None" : ~(self.left | self.right)}
        counts = {'fwd' : 0, 'bkwd' : 0}
        for step in self.simulation.storage.steps:
            step.active.sanity_check()  # traj is in ensemble
            traj = step.active[0].trajectory
            traj_str = traj.summarize_by_volumes_str(state_label)
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
        assert_equal(counts['fwd'] + counts['bkwd'], 10)

    def test_forward_only_committor(self):
        raise SkipTest

    def test_backward_only_committor(self):
        raise SkipTest
