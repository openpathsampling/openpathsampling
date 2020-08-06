from nose.tools import (assert_equal, assert_not_equal, assert_almost_equal,
                        raises, assert_true)
from nose.plugins.skip import Skip, SkipTest
from .test_helpers import assert_items_equal

import openpathsampling as paths
import openpathsampling.engines as peng
from openpathsampling.engines.toy import ToySnapshot

from openpathsampling.engines.external_engine import *

import numpy as np

import psutil
import shlex

import time
import os
import glob
import linecache

import logging

from openpathsampling.engines.snapshot import SnapshotDescriptor

logging.getLogger('openpathsampling.ensemble').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.netcdfplus').setLevel(logging.CRITICAL)

engine_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                          "external_engine")


class ExampleExternalEngine(peng.ExternalEngine):
    """Trivial external engine for engine.c in the tests.
    """
    SnaphotClass = ToySnapshot
    InternalizedSnapshotClass = ToySnapshot
    def read_frame_from_file(self, filename, frame_num):
        # under most circumstances, start with linecache.checkcache and
        # setting the value of the first line
        linecache.checkcache(filename)
        first_line = frame_num + 1

        # create a snapshot out of lines starting with first_line... if
        # nothing exists, linecache returns '', so we return None.
        # Otherwise, try to make a snapshot and return "partial" if we fail
        line = linecache.getline(filename, first_line)
        if line is '':
            snap = None
        else:
            try:
                splitted = line.split()
                if len(splitted) == 2:
                    coords = float(splitted[0])
                    vels = float(splitted[1])
                else:
                    raise ValueError()  # force the raise we then ignore
                snap = ToySnapshot(coordinates=np.array([[coords]]),
                                   velocities=np.array([[vels]]))
            except ValueError:
                snap = "partial"
        return snap

    def write_frame_to_file(self, filename, snapshot, mode='a'):
        f = open(filename, mode)
        snapshot_text = "{pos} {vel}\n".format(pos=snapshot.xyz[0][0],
                                               vel=snapshot.velocities[0][0])
        f.write(snapshot_text)
        f.close()

    def set_filenames(self, number):
        self.input_file = self.name_prefix + str(number) + ".inp"
        self.output_file = self.name_prefix + str(number) + ".out"

    def engine_command(self):
        if self.engine_directory != "":
            engine_path = os.path.join(self.engine_directory, "engine")
        else:  # pragma: no cover
            engine_path = "engine"
        return (engine_path + " " + str(self.engine_sleep)
                + " " + str(self.output_file) + " " + str(self.input_file))

def setup_module():
    proc = psutil.Popen("make", cwd=engine_dir)
    proc.wait()

def teardown_module():
    # proc = psutil.Popen("make clean", cwd=engine_dir, shell=True)
    # proc.wait()
    for testfile in glob.glob("test*out") + glob.glob("test*inp"):
        os.remove(testfile)

class TestExternalEngine(object):
    def setup(self):
        self.descriptor = SnapshotDescriptor.construct(
            snapshot_class=ToySnapshot,
            snapshot_dimensions={'n_spatial': 1,
                                 'n_atoms': 1}
        )
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
        self.slow_engine = ExampleExternalEngine(slow_options,
                                                 self.descriptor,
                                                 self.template)
        self.fast_engine = ExampleExternalEngine(fast_options,
                                                 self.descriptor,
                                                 self.template)
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

class TestFilenameSetter(object):
    def test_default_setter(self):
        setter = FilenameSetter()
        assert setter() == 0
        assert setter() == 1

    def test_specific_number_setter(self):
        setter = FilenameSetter(100)
        assert setter() == 100
        assert setter() == 101


class TestRandomStringFilenames(object):
    def test_default_setter(self):
        setter = RandomStringFilenames()
        prefix = setter()
        assert len(prefix) == 8
        for char in prefix:
            assert char in RandomStringFilenames._allowed

    def test_trivial_setter(self):
        setter = RandomStringFilenames(2, 'a')
        assert setter() == 'aa'
