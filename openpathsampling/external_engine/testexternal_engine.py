
from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        assert_almost_equal, raises)
from nose.plugins.skip import Skip, SkipTest
import openpathsampling as paths

from external_engine import *

class testExternalEngine(object):
    def setUp(self):
        options = {}
        template = paths.Snapshot()
        self.eng = ExternalEngine(options, template)

    def test_current_snapshot(self):
        raise SkipTest

    def test_read_frame_from_file(self):
        raise SkipTest

    def test_start_stop(self):
        raise SkipTest
