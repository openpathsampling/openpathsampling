from openpathsampling.experimental.simstore.custom_json import JSONCodec
from openpathsampling.experimental.simstore.attribute_handlers import (
    AttributeHandler, NDArrayHandler
)
import re

try:
    import simtk.unit
except ImportError:
    HAS_SIMTK = False
else:
    HAS_SIMTK = True
    # note that we need to force simstore to not treat this as an iterable
    from openpathsampling.experimental.simstore.class_lookup import \
        is_storage_iterable
    is_storage_iterable.force_false(simtk.unit.Quantity)

### JSON SERIALIZATION ###################################################

def unit_to_dict(obj):
    return {p.name: int(power)
            for p, power in obj.iter_base_or_scaled_units()}

def unit_from_dict(dct):
    unit = simtk.unit.Unit({})
    for u_name, u_power in dct.items():
        unit *= getattr(simtk.unit, u_name) ** u_power
    return unit

def quantity_to_dict(obj):
    return {'value': obj.value_in_unit(obj.unit),
            '__simtk_unit__': unit_to_dict(obj.unit)}

def quantity_from_dict(dct):
    unit = unit_from_dict(dct['__simtk_unit__'])
    return dct['value'] * unit


if HAS_SIMTK:
    simtk_quantity_codec = JSONCodec(
        cls=simtk.unit.Quantity,
        to_dict=quantity_to_dict,
        from_dict=quantity_from_dict,
        is_my_dict=lambda x: '__simtk_unit__' in x
    )

### DIRECT SERIALIZATION #################################################
_simtk_re = re.compile("simtk\((.*)\)\*((ndarray|float).*)")


def simtk_unit_from_string(unit_str):
    # TODO: add safety checks; parse the AST and ensure that all attributes
    # are of `unit` and that all operations are mul/div/pow
    from simtk import unit
    return eval(unit_str, {'unit': unit})


class SimtkQuantityHandler(AttributeHandler):
    # NOTE: only supports ndarrays and floats for now
    def __init__(self, type_info):
        super().__init__(type_info)
        self.unit_str, self.wrapped_type = type_info
        self.unit = simtk_unit_from_string(self.unit_str)
        if self.wrapped_type == 'float':
            self.inner_serialize = lambda x: x
            self.inner_deserialize = lambda x, _: x
            self.backend_type = 'float'
            self.type_size = None
        else:
            np_handler = NDArrayHandler.from_type_string(self.wrapped_type)
            self.inner_serialize = np_handler.serialize
            self.inner_deserialize = np_handler.deserialize
            self.backend_type = np_handler.backend_type
            self.type_size = np_handler.type_size

    @staticmethod
    def is_my_type(type_str):
        m_simtk_quantity = _simtk_re.match(type_str)
        if m_simtk_quantity:
            unit = m_simtk_quantity.group(1)
            wrapped_type = m_simtk_quantity.group(2)
            return unit, wrapped_type

    def serialize(self, obj):
        unwrapped = obj.value_in_unit(self.unit)
        return self.inner_serialize(unwrapped)

    def deserialize(self, data, caches=None):
        unwrapped = self.inner_deserialize(data, caches)
        return self.unit * unwrapped
