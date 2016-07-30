from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        assert_almost_equal, raises)
from nose.plugins.skip import Skip, SkipTest
from test_helpers import (
    true_func, assert_equal_array_array, make_1d_traj, data_filename
)

import openpathsampling as paths

import logging
logging.getLogger('openpathsampling.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.ensemble').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.netcdfplus').setLevel(logging.CRITICAL)

class testInterfaceSet(object):
    def setup(self):
        pass

    def test_get_lambda(self):
        raise SkipTest

class testGenericVolumeInterfaceSet(object):
    def setup(self):
        pass

    def test_determine_direction(self):
        raise SkipTest

    def test_prep_minvals_maxvals(self):
        raise SkipTest

    def test_bad_minvals_maxvals(self):
        raise SkipTest

class testVolumeInterfaceSet(object):
    def setup(self):
        pass

    def test_get_lambda(self):
        raise SkipTest

    def test_len(self):
        raise SkipTest

    def test_getitem(self):
        raise SkipTest

    def test_iter(self):
        raise SkipTest

    def test_contains(self):
        raise SkipTest

    def test_reversed(self):
        raise SkipTest

    def test_new_interface(self):
        raise SkipTest
    

class testPeriodicVolumeInterfaceSet(object):
    def setup(self):
        pass

    def test_get_lambda(self):
        raise SkipTest

    def test_len(self):
        raise SkipTest

    def test_getitem(self):
        raise SkipTest

    def test_iter(self):
        raise SkipTest

    def test_contains(self):
        raise SkipTest

    def test_reversed(self):
        raise SkipTest

    def test_new_interface(self):
        raise SkipTest
    
