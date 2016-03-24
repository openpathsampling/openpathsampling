from test_helpers import raises_with_message_like
from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        raises, assert_almost_equal, assert_true)
from nose.plugins.skip import SkipTest

from openpathsampling.pathsimulator import *
import openpathsampling as paths
import openpathsampling.engines.toy as toys
import numpy as np

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
        left = paths.CVRangeVolume(cv, float("-inf"), -1.0)
        right = paths.CVRangeVolume(cv, 1.0, float("inf"))

        randomizer = paths.RandomVelocities(beta=1.0)

        self.simulation = CommittorSimulation(storage=None,
                                              engine=engine,
                                              states=[left, right],
                                              randomizer=randomizer,
                                              initial_snapshots=snap0)
        

    def test_initialization(self):
        sim = self.simulation  # convenience
        assert_equal(len(sim.initial_snapshots), 1)
        assert_true(isinstance(sim.mover, paths.RandomChoiceMover))
        raise SkipTest

    def test_committor_run(self):
        raise SkipTest

    def test_forward_only_committor(self):
        raise SkipTest

    def test_backward_only_committor(self):
        raise SkipTest
