
from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        assert_almost_equal, raises)
from nose.plugins.skip import Skip, SkipTest
import openpathsampling as paths

from external_engine import *

import psutil

import time
import os

def setUp():
    # TODO: run Makefile
    pass

class testExternalEngine(object):
    def setUp(self):
        slow_options = {
            'n_frames_max' : 10000, 
            'engine_sleep' : 100,
            'name_prefix' : "test"
        }
        fast_options = {
            'n_frames_max' : 10000, 
            'engine_sleep' : 0,
            'name_prefix' : "test"
        }
        self.template = paths.Snapshot()
        self.slow_engine = ExternalEngine(slow_options, self.template)
        self.fast_engine = ExternalEngine(fast_options, self.template)
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
        assert_equal(eng.proc.status(), 'running') # zombies also run

        # stop it; check that it isn't running
        eng.stop(self.template)
        assert_equal(eng.proc.is_running(), False)

    def test_read_frame_from_file(self):
        eng = self.slow_engine
        testf = open('testf1.data', 'w')
        testf.write("1.0\n2.0\n3.0\n")
        testf.close()
        snap2 = eng.read_frame_from_file("testf1.data", 2)
        assert_equal(snap2.xyz[0][0], 2.0)
        snap1 = eng.read_frame_from_file("testf1.data", 1)
        assert_equal(snap1.xyz[0][0], 1.0)
        snap3 = eng.read_frame_from_file("testf1.data", 3)
        assert_equal(snap3.xyz[0][0], 3.0)
        snap4 = eng.read_frame_from_file("testf1.data", 4)
        assert_equal(snap4, None)
        os.remove('testf1.data')

    def test_read_frame_while_writing_file(self):
        eng = self.slow_engine
        testf = open('testf2.data', 'w')
        testf.write("6.0\ninvalid")
        testf.close()
        snap1 = eng.read_frame_from_file("testf2.data", 1)
        assert_equal(snap1.xyz[0][0], 6.0)
        snap2 = eng.read_frame_from_file("testf2.data", 2)
        assert_equal(snap2, "partial")
        os.remove('testf2.data')

    def test_slow_run(self):
        # generate traj in LengthEnsemble if frames only come every 100ms
        raise SkipTest

    def test_fast_run(self):
        # generate traj in LengthEnsemble if frames come as fast as possible
        raise SkipTest


