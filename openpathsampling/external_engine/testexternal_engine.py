
from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        assert_almost_equal, raises)
from nose.plugins.skip import Skip, SkipTest
import openpathsampling as paths

from external_engine import *

class testExternalEngine(object):
    def setUp(self):
        slow_options = {
            'n_frames_max' : 10000, 
            'engine_command' : "engine 100 engine.out",
            'trajectory_file' : "engine.out"
        }
        fast_options = {
            'n_frames_max' : 10000, 
            'engine_command' : "engine 0 engine.out",
            'trajectory_file' : "engine.out"
        }
        self.template = paths.Snapshot()
        self.slow_engine = ExternalEngine(slow_options, self.template)
        self.fast_engine = ExternalEngine(fast_options, self.template)
        self.ensemble = paths.LengthEnsemble(5)

    def test_start_stop(self):
        self.slow_engine.start(self.template)
        # check that it is running
        self.slow_engine.stop(self.template)
        # check that it has stopped
        raise SkipTest

    def test_read_frame_from_file(self):
        raise SkipTest

    def test_slow_run(self):
        # generate traj in LengthEnsemble if frames only come every 100ms
        raise SkipTest

    def test_fast_run(self):
        # generate traj in LengthEnsemble if frames come as fast as possible
        raise SkipTest

