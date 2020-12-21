import pytest
from .type_ident import *

from .test_utils import MockUUIDObject

_TYPES = ['int', 'float', 'str', 'ndarray', 'bool', 'uuid']

class TestStandardTyping(object):
    # this tests everything in STANDARD_TYPING; it's a little more
    # integration test than unit, but the individual tests ensure that each
    # unit is tested
    def setup(self):
        self.objects = {
            'int': ('int', 5),
            'float': ('float', 2.3),
            'str': ('str', "foo"),
            'ndarray': ('ndarray.float64(2,3)',
                        np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])),
            'bool': ('bool', True),
            'uuid': ('uuid', MockUUIDObject(name='int', normal_attr=5)),
        }

    @pytest.mark.parametrize('example', _TYPES)
    def test_identify(self, example):
        string, obj = self.objects[example]
        assert STANDARD_TYPING.identify(obj) == string

    @pytest.mark.parametrize('example', _TYPES)
    def test_parsing(self, example):
        string, _ = self.objects[example]
        assert STANDARD_TYPING.parse(string) == example

    def test_error_register_existing(self):
        with pytest.raises(RuntimeError):
            STANDARD_TYPING.register(int_type_id)

    def test_error_identify(self):
        with pytest.raises(TypeIdentificationError):
            STANDARD_TYPING.identify(object())

    def test_error_parse(self):
        with pytest.raises(TypeStringError):
            STANDARD_TYPING.parse('foo')
