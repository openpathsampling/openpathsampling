from nose.tools import assert_equal, assert_not_equal, raises
from nose.plugins.skip import SkipTest
from test_helpers import CallIdentity, prepend_exception_message

from opentis.globalstate import *

class testGlobalState(object):
    def setup(self):
        pass

