from __future__ import absolute_import
from nose.tools import (assert_equal, assert_not_equal, assert_almost_equal,
                        raises)
from nose.plugins.skip import Skip, SkipTest

from ..test_helpers import (
    true_func, assert_equal_array_array, make_1d_traj, data_filename,
    assert_items_equal
)

import openpathsampling as paths

from openpathsampling.pathsimulators.hooks import *

import logging
logging.getLogger('openpathsampling.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.ensemble').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.netcdfplus').setLevel(logging.CRITICAL)

class TestPathSimulatorHook(object):
    def setup(self):
        pass

    def test_attach_hooks(self):
        pass

class TestStorageHook(object):
    def setup(self):
        pass

    def test_before_simulation(self):
        pass

    def test_after_step(self):
        pass

    def test_after_simulation(self):
        pass

class TestShootFromSnapshotsOutputHook(object):
    def setup(self):
        pass

    def test_before_simulation(self):
        pass

    def test_before_step(self):
        pass
