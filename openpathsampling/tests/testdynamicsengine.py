import openpathsampling as paths

from nose.tools import (assert_equal, assert_not_equal, raises)
from nose.plugins.skip import SkipTest
from test_helpers import make_1d_traj, raises_with_message_like

class StupidEngine(paths.engines.DynamicsEngine):
    _default_options = {'random_option': False}

    @property
    def bad_property(self):
        obj = object()
        return obj.b

    @property
    def property_recovers(self):
        if not hasattr(self, 'attempted'):
            self.attempted = True
            raise AttributeError("Internal error")
        return self.attempted

class testDynamicsEngine(object):
    def setup(self):
        options = {'n_frames_max' : 100, 'random_option' : True}
        template = make_1d_traj([1.0])[0]
        self.engine = paths.engines.DynamicsEngine(options, template)
        self.stupid = StupidEngine(options, template)

    def test_getattr_from_options(self):
        assert_equal(self.stupid.random_option, True)

    @raises_with_message_like(AttributeError, 
                              "'object' object has no attribute 'b'")
    def test_getattr_property_fails(self):
        self.stupid.bad_property

    @raises_with_message_like(AttributeError, 
                              "Unknown problem occurred in property")
    def test_getattr_property_recovers(self):
        self.stupid.property_recovers

    @raises_with_message_like(AttributeError,
                              "'StupidEngine' has no attribute 'foo'" +  
                              ", nor does its options dictionary")
    def test_getattr_does_not_exist(self):
        self.stupid.foo
