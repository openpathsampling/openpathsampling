
from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        assert_almost_equal, raises)
from nose.plugins.skip import Skip, SkipTest
import openpathsampling as paths

from external_engine import *

import psutil

import time

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
        #print get_all_pids()
        # check that it isn't running yet
        try:
            assert_equal(eng.proc.is_running(), False)
        except AttributeError:
            pass # if eng.proc doesn't exist, then it isn't running

        # start it; check that it is running
        eng.start(self.template)
        assert_equal(eng.proc.is_running(), True)
        assert_equal(eng.proc.status(), 'running')

        # stop it; check that it isn't running
        eng.stop(self.template)
        assert_equal(eng.proc.is_running(), False)

    def test_read_frame_from_file(self):
        raise SkipTest

    def test_slow_run(self):
        # generate traj in LengthEnsemble if frames only come every 100ms
        raise SkipTest

    def test_fast_run(self):
        # generate traj in LengthEnsemble if frames come as fast as possible
        raise SkipTest


