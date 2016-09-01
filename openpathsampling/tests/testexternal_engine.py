from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        assert_almost_equal, raises, assert_true)
from nose.plugins.skip import Skip, SkipTest

import openpathsampling as paths
import openpathsampling.engines as peng

import numpy as np

import psutil
import shlex

import time
import os
import glob

import logging

logging.getLogger('openpathsampling.ensemble').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.netcdfplus').setLevel(logging.CRITICAL)

engine_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                          "external_engine")

def setUp():
    proc = psutil.Popen("make", cwd=engine_dir)
    proc.wait()

def teardown():
    # proc = psutil.Popen("make clean", cwd=engine_dir, shell=True)
    # proc.wait()
    for testfile in glob.glob("test*out") + glob.glob("test*inp"):
        os.remove(testfile)

class testExternalEngine(object):
    def setUp(self):
        slow_options = {
            'n_frames_max' : 10000, 
            'engine_sleep' : 100,
            'name_prefix' : "test",
            'engine_directory' : engine_dir
        }
        fast_options = {
            'n_frames_max' : 10000, 
            'engine_sleep' : 0,
            'name_prefix' : "test",
            'engine_directory' : engine_dir
        }
        self.template = peng.toy.Snapshot(coordinates=np.array([[0.0]]),
                                          velocities=np.array([[1.0]]))
        self.slow_engine = peng.ExternalEngine(slow_options, self.template)
        self.fast_engine = peng.ExternalEngine(fast_options, self.template)
        self.ensemble = paths.LengthEnsemble(5)

    def test_start_stop(self):
        eng = self.fast_engine
        # check that it isn't running yet
        try:
            assert_equal(eng.proc.is_running(), False)
        except AttributeError:
            pass # if eng.proc doesn't exist, then it isn't running

        # start it; check that it is running
        eng.start(self.template)
        assert_equal(eng.proc.is_running(), True)
        # zombies also run
        assert_not_equal(eng.proc.status(), psutil.STATUS_ZOMBIE) 

        # stop it; check that it isn't running
        eng.stop(None)
        assert_equal(eng.proc.is_running(), False)

    def test_read_frame_from_file(self):
        eng = self.slow_engine
        testf = open('testf1.data', 'w')
        testf.write("1.0 1.0\n2.0 1.0\n3.0 1.0\n")
        testf.close()
        snap2 = eng.read_frame_from_file("testf1.data", 1)
        assert_equal(snap2.xyz[0][0], 2.0)
        snap1 = eng.read_frame_from_file("testf1.data", 0)
        assert_equal(snap1.xyz[0][0], 1.0)
        snap3 = eng.read_frame_from_file("testf1.data", 2)
        assert_equal(snap3.xyz[0][0], 3.0)
        snap4 = eng.read_frame_from_file("testf1.data", 3)
        assert_equal(snap4, None)
        os.remove('testf1.data')

    def test_read_frame_while_writing_file(self):
        eng = self.slow_engine
        testf = open('testf2.data', 'w')
        testf.write("6.0 1.0\ninvalid")
        testf.close()
        snap1 = eng.read_frame_from_file("testf2.data", 0)
        assert_equal(snap1.xyz[0][0], 6.0)
        snap2 = eng.read_frame_from_file("testf2.data", 1)
        assert_equal(snap2, "partial")
        os.remove('testf2.data')

        testf = open('testf3.data', 'w')
        testf.write("6.0 ")
        testf.close()
        snap3 = eng.read_frame_from_file("testf3.data", 0)
        assert_equal(snap3, "partial")
        os.remove('testf3.data')

    def test_generate_next_frame(self):
        eng = self.slow_engine
        eng.start(self.template)
        traj = eng.generate_next_frame()
        eng.stop(traj)

    def test_slow_run(self):
        # generate traj in LengthEnsemble if frames only come every 100ms
        self.slow_engine.initialized = True
        traj = self.slow_engine.generate(self.template, 
                                         [self.ensemble.can_append])
        assert_equal(len(traj), 5)

    def test_fast_run(self):
        # generate traj in LengthEnsemble if frames come as fast as possible
        self.fast_engine.initialized = True
        traj = self.fast_engine.generate(self.template, 
                                         [self.ensemble.can_append])
        assert_equal(len(traj), 5)

    def test_in_shooting_move(self):
        for testfile in glob.glob("test*out") + glob.glob("test*inp"):
            os.remove(testfile)
        ens10 = paths.LengthEnsemble(10)
        init_traj = self.fast_engine.generate(self.template,
                                              [ens10.can_append])
        assert_equal(ens10(init_traj), True)
        init_conds = paths.SampleSet([
            paths.Sample(replica=0, ensemble=ens10, trajectory=init_traj)
        ])
        shooter = paths.OneWayShootingMover(ensemble=ens10,
                                            selector=paths.UniformSelector(),
                                            engine=self.fast_engine)
        prev_sample_set = init_conds
        default_traj = [[[0.0]], [[1.0]], [[2.0]], [[3.0]], [[4.0]],
                        [[5.0]], [[6.0]], [[7.0]], [[8.0]], [[9.0]]]
        assert_items_equal(init_conds[0].trajectory.xyz, default_traj)
        for step in range(10):
            assert_equal(len(prev_sample_set), 1)
            change = shooter.move(prev_sample_set)
            new_sample_set = prev_sample_set.apply_samples(change.results)
            assert_items_equal(new_sample_set[0].trajectory.xyz,
                               default_traj)
            prev_traj = prev_sample_set[0].trajectory
            new_traj = new_sample_set[0].trajectory
            shared = prev_traj.shared_configurations(new_traj)
            assert_true(0 < len(list(shared)) < len(new_traj))
            prev_sample_set = new_sample_set

        for testfile in glob.glob("test*out") + glob.glob("test*inp"):
            os.remove(testfile)

