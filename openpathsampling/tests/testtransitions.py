from nose.tools import assert_equal, assert_not_equal, raises
from nose.plugins.skip import SkipTest
from test_helpers import CallIdentity, prepend_exception_message

import openpathsampling as paths
from openpathsampling.ensemble import *

import logging
logging.getLogger('opentis.analysis.tis_analysis').setLevel(logging.DEBUG)
logging.getLogger('opentis.initialization').setLevel(logging.CRITICAL)

class testTISTransition(object):
    def setup(self):
        pass

    def test_initialization(self):
        pass

    def test_ensemble_statistics(self):
        pass

