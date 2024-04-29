import pytest

from openpathsampling.netcdfplus.base import *

class TestStorableNamedObject(object):
    def setup_method(self):
        self.obj = StorableNamedObject()

    def test_named(self):
        assert not self.obj.is_named
        obj = self.obj.named('foo')
        assert obj is self.obj
        assert self.obj.is_named
        assert self.obj.name == 'foo'

    def test_bad_name_type(self):
        with pytest.raises(TypeError, match="Invalid name type"):
            self.obj.named(1)
