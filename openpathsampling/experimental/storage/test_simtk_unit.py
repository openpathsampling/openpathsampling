import pytest
import numpy as np
import tempfile
import os

from ..simstore.custom_json import JSONSerializerDeserializer, DEFAULT_CODECS
from ..storage import Storage

from .simtk_unit import *

try:
    from simtk import unit
except ImportError:
    HAS_SIMTK = False
else:
    HAS_SIMTK = True

class TestSimtkUnitCodec(object):
    def setup(self):
        pytest.importorskip('simtk.unit')
        my_unit = unit.nanometer / unit.picosecond**2
        self.values = {
            'float': 1.0 * my_unit,
            'array': np.array([1.0, 2.0]) * my_unit,
        }
        self.serialization = JSONSerializerDeserializer(
            DEFAULT_CODECS + [simtk_quantity_codec]
        )

    @pytest.mark.parametrize('obj_type', ['float', 'array'])
    def test_serialization_cycle(self, obj_type):
        obj = self.values[obj_type]
        ser = self.serialization.serializer(obj)
        deser = self.serialization.deserializer(ser)
        reser = self.serialization.serializer(deser)
        if obj_type == 'array':
            np.testing.assert_array_equal(obj, deser)
        else:
            assert obj == deser
        assert ser == reser


class TestSimtkQuantityHandler(object):
    def setup(self):
        pytest.importorskip('simtk.unit')
        self.handlers = {
            'float': SimtkQuantityHandler(
                ('unit.nanometer/unit.picosecond**2', 'float')
            ),
            'array': SimtkQuantityHandler(
                ('unit.nanometer', 'ndarray.float32(2,3)')
            ),
        }
        self.objects = {
            'float': 1.0 * unit.nanometer / unit.picosecond**2,
            'array': np.array([[1.0, 2.0, 3.0],
                               [4.0, 5.0, 6.0]]) * unit.nanometer,
        }

    @pytest.mark.parametrize('type_str, expected', [
        (
            'simtk(unit.nanometer/unit.picosecond**2)*float',
             ('unit.nanometer/unit.picosecond**2', 'float')
        ), (
            'simtk(unit.nanometer)*ndarray.float32(3,3)',
            ('unit.nanometer', 'ndarray.float32(3,3)')
        ),
    ])
    def test_is_my_type(self, type_str, expected):
        assert SimtkQuantityHandler.is_my_type(type_str) == expected

    @pytest.mark.parametrize('obj_type', ['float', 'array'])
    def test_serialization_cycle(self, obj_type):
        handler = self.handlers[obj_type]
        obj = self.objects[obj_type]
        ser = handler.serialize(obj)
        deser = handler.deserialize(ser)
        reser = handler.serialize(deser)

        assert ser == reser
        if obj_type == 'array':
            np.testing.assert_array_equal(obj, deser)
        else:
            assert obj == deser
        assert obj.unit == deser.unit


@pytest.mark.parametrize('obj_type', ['float', 'array'])
def test_tag_simtk(obj_type):
    # this is mainly a smoke test against regression that we couldn't tag
    # simtk quantities
    pytest.importorskip('simtk.unit')
    value = {'float': 2.0, 'array': np.array([1.0, 2.0])}[obj_type]
    quantity = value * unit.nanometer
    with tempfile.TemporaryDirectory() as tmpdir:
        filename = os.path.join(tmpdir, "test.db")
        storage = Storage(filename, mode='w')
        storage.tags['foo'] = quantity
        storage2 = Storage(filename, mode='r')
        reloaded = storage2.tags['foo']
        assert isinstance(reloaded, unit.Quantity)
        assert reloaded is not quantity
        if obj_type == 'array':
            np.testing.assert_array_equal(quantity, reloaded)
        else:
            assert quantity == reloaded
        assert quantity.unit == reloaded.unit

