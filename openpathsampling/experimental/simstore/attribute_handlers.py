from .my_types import parse_ndarray_type
import numpy as np


class AttributeHandler(object):
    def __init__(self, type_info):
        self.type_info = type_info

    @staticmethod
    def is_my_type(type_str):
        """Returns type info (possibly tuple) if true, None if false"""
        raise NotImplementedError()

    @classmethod
    def from_type_string(cls, type_str):
        type_info = cls.is_my_type(type_str)
        if type_info:
            return cls(type_info)

    def serialize(self, obj):
        return obj

    def deserialize(self, data, caches=None):
        return data


class StandardHandler(AttributeHandler):
    standard_types = ['str', 'int', 'float', 'function']
    def __init__(self, type_info):
        super().__init__(type_info)
        self.backend_type = type_info
        self.type_size = None

    @classmethod
    def is_my_type(cls, type_str):
        if type_str in cls.standard_types:
            return type_str


# TODO: these have not yet been implemented, but this is how it should all
# work -- everything should be managed through attribute handlers, instead
# of separate functions for serialize/deserialize
class UUIDHandler(StandardHandler):
    standard_types = ['uuid', 'lazy']
    def serialize(self, obj):
        pass

    def deserialize(self, data, caches=None):
        pass


class ListUUIDHandler(StandardHandler):
    standard_types = ['list_uuid']
    def serialize(self, obj):
        pass

    def deserialize(self, data, caches=None):
        pass


class JSONObjHandler(StandardHandler):
    standard_types = ['json_obj']
    def serialize(self, obj):
        pass

    def deserialize(self, data, caches=None):
        pass


class NDArrayHandler(AttributeHandler):
    def __init__(self, type_info):
        super().__init__(type_info)
        self.dtype, self.shape = type_info
        self.backend_type = 'ndarray'
        self.type_size = None  # TODO: change this based on dtype/shape

    @classmethod
    def is_my_type(cls, type_str):
        return parse_ndarray_type(type_str)

    def serialize(self, obj):
        return obj.astype(dtype=self.dtype, copy=False).tostring()

    def deserialize(self, data, caches=None):
        return np.fromstring(data, dtype=self.dtype).reshape(self.shape)

DEFAULT_HANDLERS = [NDArrayHandler, StandardHandler]
