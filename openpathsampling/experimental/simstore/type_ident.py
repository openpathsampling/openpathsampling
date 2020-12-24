import re
import numbers
import numpy as np
from .serialization_helpers import has_uuid

class TypeStringError(RuntimeError):
    pass

class TypeIdentificationError(RuntimeError):
    pass

class TypingManager(object):
    def __init__(self, type_ids):
        self._type_ids = {}
        for type_id in type_ids:
            self.register(type_id)

    def register(self, type_id, force=False):
        if not force and type_id.name in self._type_ids:
            raise RuntimeError("Type identifier for %s already exists. "
                               "Use `force=True` to override."
                               % type_id.name)
        self._type_ids[type_id.name] = type_id

    def _typing_func(self, inp, func_name):
        for type_id in self._type_ids.values():
            func = getattr(type_id, func_name)
            result = func(inp)
            if result is not None:
                break
        return result

    def parse(self, string):
        result = self._typing_func(string, 'parse')
        if result is None:
            raise TypeStringError("Unable to parse type string: ", string)
        return result

    def identify(self, obj):
        result = self._typing_func(obj, 'identify')
        if result is None:
            raise TypeIdentificationError("Unable to identify backend type "
                                          "for object of type "
                                          + str(type(obj)))
        return result


class StandardTypeIdentifier(object):
    """TypeIdentifier for most simple types.

    Simple types are types where there is only one string identifier for the
    type, and where that string can be determined using ``isinstance``.
    """
    # NOTE: Current implementation uses isinstance, which might be slow.
    # Performance could be improved by using ClassIsSomething, if needed.
    # However, expectation is for identify to be rarely used.
    def __init__(self, name, cls):
        self.name = name
        self.cls = cls

    def parse(self, string):
        if string == self.name:
            return self.name

    def identify(self, obj):
        if isinstance(obj, self.cls):
            return self.name


class NumpyTypeIdentifier(object):
    ndarray_re = re.compile(
        "ndarray\.(?P<dtype>[a-z0-9]+)(?P<shape>\([0-9\,\ ]+\))"
    )
    def __init__(self):
        self.name = 'ndarray'

    def parse(self, string):
        m_ndarray = self.ndarray_re.match(string)
        if m_ndarray:
            return self.name

    @staticmethod
    def identify(obj):
        if isinstance(obj, np.ndarray):
            dtype = str(obj.dtype)
            shape = "(" + ",".join(str(a) for a in obj.shape) + ")"
            return "ndarray." + dtype + shape


class UUIDTypeIdentifier(StandardTypeIdentifier):
    def identify(self, obj):
        if has_uuid(obj):
            return self.name


int_type_id = StandardTypeIdentifier('int', numbers.Integral)
float_type_id = StandardTypeIdentifier('float', numbers.Real)
str_type_id = StandardTypeIdentifier('str', str)
bool_type_id = StandardTypeIdentifier('bool', bool)
numpy_type_id = NumpyTypeIdentifier()
uuid_type_id = UUIDTypeIdentifier('uuid', cls=None)


STANDARD_TYPING = TypingManager([bool_type_id, int_type_id, float_type_id,
                                 str_type_id, numpy_type_id, uuid_type_id])
