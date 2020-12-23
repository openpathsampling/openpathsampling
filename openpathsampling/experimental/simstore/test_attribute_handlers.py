import pytest

import numpy as np

from .attribute_handlers import *

class TestStandardHandler(object):
    def setup(self):
        self.obj = {'str': 'foo',
                    'int': 42,
                    'float': 3.14159,
                    'function': lambda x: x}

    @pytest.mark.parametrize('type_str, expected', [
        ('str', 'str'), ('int', 'int'), ('float', 'float'),
        ('function', 'function'), ('ndarray.float32(3,2)', None),
    ])
    def test_is_my_type(self, type_str, expected):
        assert StandardHandler.is_my_type(type_str) == expected

    @pytest.mark.parametrize('type_str', ['str', 'int', 'float', 'function'])
    def test_serialization_cycle(self, type_str):
        handler = StandardHandler(type_str)
        data = self.obj[type_str]
        ser = handler.serialize(data)
        deser = handler.deserialize(ser)
        reser = handler.serialize(deser)
        assert data == deser
        assert ser == reser


class TestNDArrayHandler(object):
    def setup(self):
        self.data = np.array([[0.0, 1.0, 2.0], [3.0, 4.0, 5.0]])
        self.ndarray_handler = NDArrayHandler(('float32', (2,3)))

    @pytest.mark.parametrize('type_str,expected', [
        ('str', None),
        ('ndarray.float32(3, 2)', (np.float32, (3, 2))),
        ('ndarray.float32(3,2)', (np.float32, (3, 2))),
    ])
    def test_is_my_type(self, type_str, expected):
        assert NDArrayHandler.is_my_type(type_str) == expected

    def test_serialization_cycle(self):
        ser = self.ndarray_handler.serialize(self.data)
        deser = self.ndarray_handler.deserialize(ser)
        np.testing.assert_equal(deser, self.data)
        reser = self.ndarray_handler.serialize(deser)
        assert ser == reser
