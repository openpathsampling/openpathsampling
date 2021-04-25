from .my_types import parse_ndarray_type
import numpy as np


class AttributeHandler(object):
    """Abstract object to handler a given attribute type.

    Each attribute in a schema entry will generate an instance of this for
    the appropriate type. Different subclasses handle different types that
    are supported by SimStore.

    Note that you'll usually use this with the :method:`.from_type_string`
    constructor, which returns None if the type string doesn't match what
    this class handles.
    """
    def __init__(self, type_info):
        self.type_info = type_info

    @staticmethod
    def is_my_type(type_str):
        """Returns type info (possibly tuple) if true, None if false"""
        raise NotImplementedError()

    @classmethod
    def from_type_string(cls, type_str):
        """Generate attribute handler based on a given type string.

        Parameters
        ----------
        type_str: str
            the specific type string in a format understood by SimStore

        Returns
        -------
        :class:`.AttributeHandler` or None:
            an attribute handler for the given type string, or None if this
            class doesn't know how to handle that type
        """
        type_info = cls.is_my_type(type_str)
        if type_info:
            return cls(type_info)

    def serialize(self, obj):
        """Serialize the object.

        Default serialization just returns the object.

        Parameters
        ----------
        obj : Any
            object to be serialized

        Returns
        -------
        Any:
            version ready to be written to disk; typically str or bytes but
            standard types are also allowed (bool, float, int, etc.)
        """
        return obj

    def deserialize(self, data, caches=None):
        """Deserialize the serialized form.

        Parameters
        ----------
        data: Any
            bytes to deserialize (usually as str or bytes)
        caches: Dict
            mapping of identifiers to known objects that might be components
            of this object; mainly used when deserializing UUID objects.
        """
        return data


class StandardHandler(AttributeHandler):
    """Attribute handler for common standard types
    """
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
    """Attribute handler for NumPy ndarrays.

    Parameters
    ----------
    type_info: Tuple
        dtype and shape of the array
    """
    def __init__(self, type_info):
        super().__init__(type_info)
        self.dtype, self.shape = type_info
        self.backend_type = 'ndarray'
        self.type_size = None  # TODO: change this based on dtype/shape

    @classmethod
    def is_my_type(cls, type_str):
        return parse_ndarray_type(type_str)

    def serialize(self, obj):
        return obj.astype(dtype=self.dtype, copy=False).tobytes()

    def deserialize(self, data, caches=None):
        return np.frombuffer(data, dtype=self.dtype).reshape(self.shape)

DEFAULT_HANDLERS = [NDArrayHandler, StandardHandler]
