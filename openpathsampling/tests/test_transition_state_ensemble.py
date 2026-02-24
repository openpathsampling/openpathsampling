import os
import numpy as np
import openpathsampling as paths
import openpathsampling.engines.toy as toys
from openpathsampling.tests.test_helpers import data_filename
from nose.tools import assert_equal, assert_is_none, assert_list_equal
from nose.tools import raises
from openpathsampling.analysis.transition_state_ensemble import *

import logging
logging.getLogger('openpathsampling.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.netcdfplus').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.ensemble').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.engines').setLevel(logging.CRITICAL)


class TestTransitionStateEnsemble(object):

    def setup(self):
        pes = toys.LinearSlope(m=[0.0], c=[0.0])  # flat line
        topology = toys.Topology(n_spatial=1, masses=[1.0], pes=pes)
        integrator = toys.LeapfrogVerletIntegrator(0.1)
        options0 = {
            'integ': integrator,
            'n_frames_max': 100000,
            'n_steps_per_frame': 5
        }
        options1 = {
            'integ': integrator,
            'n_frames_max': 5,
            'n_steps_per_frame': 1
        }

        self.engine0 = toys.Engine(options=options0, topology=topology)
        self.engine1 = toys.Engine(options=options1, topology=topology)
        self.snap0 = toys.Snapshot(coordinates=np.array([[0.0]]),
                                   velocities=np.array([[1.0]]),
                                   engine=self.engine0)
        self.snap1 = toys.Snapshot(coordinates=np.array([[0.0]]),
                                   velocities=np.array([[1.0]]),
                                   engine=self.engine1)
        cv = paths.FunctionCV("Id", lambda snap: snap.coordinates[0][0])
        self.left = paths.CVDefinedVolume(cv, float("-inf"), -1.0)
        self.right = paths.CVDefinedVolume(cv, 1.0, float("inf"))
        self.state_labels = {"Left": self.left,
                             "Right": self.right,
                             "None": ~(self.left | self.right)}

        self.randomizer = paths.NoModification()

        self.filename = data_filename("committor_test.nc")
        self.storage = paths.Storage(self.filename,
                                     mode="w")
        self.storage.save(self.snap0)

    def teardown(self):
        if os.path.isfile(self.filename):
            os.remove(self.filename)

    @raises(ValueError)
    def test_pB_sanity(self):
        TransitionStateEnsemble(None, None, None, None, None, None,
                                pB_min=1.0,
                                pB_max=0.0)

    def test_no_list_trajectories(self):
        tse = TransitionStateEnsemble([self.snap0]*3,
                                      stateA=self.left, stateB=self.right,
                                      engine=self.engine0,
                                      storage=self.storage,
                                      randomizer=self.randomizer,
                                      pB_min=0.0, pB_max=1.0,
                                      n_per_snapshot=1)
        assert_list_equal([[self.snap0]*3], tse.trajectories)

    def test_working(self):
        tse = TransitionStateEnsemble([[self.snap0]*3, [self.snap1]*4],
                                      stateA=self.left, stateB=self.right,
                                      engine=self.engine0,
                                      storage=self.storage,
                                      randomizer=self.randomizer,
                                      pB_min=0.0, pB_max=1.0,
                                      n_per_snapshot=1)

        result = tse.run()
        assert_equal(len(result.keys()), 2)

    @raises(RuntimeError)
    def test_n_frames_max_error(self):
        tse = TransitionStateEnsemble([[self.snap1]*3], stateA=self.left,
                                      stateB=self.right, engine=self.engine1,
                                      storage=self.storage,
                                      randomizer=self.randomizer,
                                      pB_min=0.0, pB_max=1.0,
                                      n_per_snapshot=1)

        tse.run()


class TestNextFrame:
    def setup(self):
        self.tse = TransitionStateEnsemble(None, None, None, None, None, None)
        self.tse.snap_frame = 5

    def test_shift_right_even(self):
        self.tse.next_frame(pB=0.0, pB_min=0.1, pB_max=1.0,
                            snap_min=0,
                            snap_frame=5,
                            snap_max=10)
        assert_equal(self.tse.snap_frame, 7)

    def test_shift_right_uneven(self):
        self.tse.next_frame(pB=0.0, pB_min=0.1, pB_max=1.0,
                            snap_min=0,
                            snap_frame=5,
                            snap_max=11
                            )
        assert_equal(self.tse.snap_frame, 8)

    def test_shift_left_even(self):
        self.tse.next_frame(pB=0.3, pB_min=0.1, pB_max=0.2,
                            snap_min=0,
                            snap_frame=5,
                            snap_max=10)
        assert_equal(self.tse.snap_frame, 2)

    def test_shift_left_uneven(self):
        self.tse.next_frame(pB=0.3, pB_min=0.1, pB_max=0.2,
                            snap_min=0,
                            snap_frame=6,
                            snap_max=10
                            )
        assert_equal(self.tse.snap_frame, 3)

    @raises(ValueError)
    def test_pB_sanity(self):
        self.tse.next_frame(pB=0.2, pB_min=0.3, pB_max=0.1,
                            snap_min=0,
                            snap_frame=5,
                            snap_max=10
                            )

    @raises(RuntimeError)
    def test_selection_left_sanity(self):
        self.tse.next_frame(pB=0.3, pB_min=0.1, pB_max=0.2,
                            snap_min=1,
                            snap_frame=2,
                            snap_max=3
                            )

    @raises(RuntimeError)
    def test_selection_right_sanity(self):
        self.tse.next_frame(pB=0.0, pB_min=0.1, pB_max=0.2,
                            snap_min=1,
                            snap_frame=2,
                            snap_max=3
                            )

    @raises(ValueError)
    def test_snap_min_max_sanity(self):
        self.tse.next_frame(pB=0.3, pB_min=0.1, pB_max=0.2,
                            snap_min=3,
                            snap_frame=2,
                            snap_max=1
                            )

    def test_accept(self):
        self.tse.next_frame(pB=0.2, pB_min=0.1, pB_max=0.3,
                            snap_min=3,
                            snap_frame=2,
                            snap_max=1
                            )
        assert_is_none(self.tse.snap_frame)

    def test_pB_min_accept(self):
        self.tse.next_frame(pB=0.1, pB_min=0.1, pB_max=0.3,
                            snap_min=3,
                            snap_frame=2,
                            snap_max=1
                            )
        assert_is_none(self.tse.snap_frame)

    def test_pB_max_accept(self):
        self.tse.next_frame(pB=0.3, pB_min=0.1, pB_max=0.3,
                            snap_min=3,
                            snap_frame=2,
                            snap_max=1
                            )
        assert_is_none(self.tse.snap_frame)
