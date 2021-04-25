from builtins import object
from openpathsampling.rng import *


class TestRandom(object):
    def test_default_rng(self):
        rng = default_rng()
        # Test if we can grab an integer (and rng excludes endpoint)
        i = rng.integers(1)
        assert i == 0

    def test_identities(self):
        rng1 = default_rng()
        rng2 = default_rng()
        # These should be te same object
        assert rng1 is rng2
