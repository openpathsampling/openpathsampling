from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        raises, assert_almost_equal)
from nose.plugins.skip import SkipTest
from test_helpers import make_1d_traj

import openpathsampling as paths

import logging
logging.getLogger('openpathsampling.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.netcdfplus').setLevel(logging.CRITICAL)

class testNoModification(object):
    def setup(self):
        pass

    # test methods from the abstract base class in here:
    def test_extract_subset(self):
        pass

    def test_apply_to_subset(self):
        pass

    def test_call(self):
        pass

class testRandomizeVelocities(object):
    def setup(self):
        pass

    def test_call(self):
        pass
