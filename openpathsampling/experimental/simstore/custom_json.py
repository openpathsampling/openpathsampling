import json
import functools
import inspect
from collections import namedtuple, defaultdict
from .tools import none_to_default
from .serialization_helpers import (
    has_uuid, replace_uuid, encode_uuid, get_uuid, set_uuid
)
from .serialization_helpers import do_import, from_dict_with_uuids

class SimulationObjectSerialization(object):
    def __init__(self, json_encoder, json_decoder):
        self.json_encoder = json_encoder
        self.json_decoder = json_decoder

    def serializer(self, obj):
        return {'uuid': get_uuid(obj),
                'json': self.json_encoder(obj)}

    def deserialize(self, uuid, table_row, cache_list):
        # TODO: is uuid input necessary here?
        dct = self.json_decoder(table_row['json'])
        cls = do_import(dct.pop('__module__'), dct.pop('__class__'))

        # TODO: this is a hack around some objects having name params
        has_name_param = 'name' in inspect.signature(cls).parameters
        name = None if has_name_param else dct.pop('name', None)
        # should just be dct.pop('name', None)

        dct = from_dict_with_uuids(dct, cache_list)
        obj = cls.from_dict(dct)
        set_uuid(obj, uuid)
        if name:
            obj.name = name
        return obj


class JSONSerializerDeserializer(object):
    """
    Tools to serialize and deserialize objects as JSON.

    This wrapper object is necessary so that we can register new codecs
    after the original initialization.

    Parameters
    ----------
    codecs : list of :class:`.JSONCodec`s
        codecs supported
    """
    def __init__(self, codecs, named_codecs=None):
        self.named_codecs = none_to_default(named_codecs, {})
        self.codecs = []
        for codec in codecs:
            self.add_codec(codec)
        self._set_serialization()

    def _set_serialization(self):
        encoder, decoder = custom_json_factory(self.codecs)
        self._serializer = functools.partial(json.dumps, cls=encoder)
        self._deserializer = functools.partial(json.loads, cls=decoder)
        self._sim_serialization = SimulationObjectSerialization(
            self._serializer, self._deserializer
        )

    def add_codec(self, codec):
        """Add a new codec to the supported codecs

        Parameters
        ----------
        codec : :class:`.JSONCodec`
            codec to add
        """
        if codec in self.codecs:
            return

        if codec is not None:
            self.codecs.append(codec)

        self._set_serialization()

    def replace_named_codec(self, codec_name, codec):
        self.named_codecs[codec_name] = codec
        self.add_codec(None)

    def serializer(self, obj):
        return self._serializer(obj)

    def deserializer(self, string):
        return self._deserializer(string)

    def simobj_serializer(self, obj):
        return self._sim_serialization.serializer(obj)

    def simobj_deserializer(self, uuid, table_row, cache_list):
        return self._sim_serialization.deserialize(uuid, table_row,
                                                   cache_list)


def custom_json_factory(coding_methods):
    """Create JSONEncoder/JSONDecoder for special types
    """
    class CustomJSONEncoder(json.JSONEncoder):
        def default(self, obj):
            for coding_method in coding_methods:
                result = coding_method.default(obj)
                if result is not obj:
                    return result
            return json.JSONEncoder.default(self, obj)

    class CustomJSONDecoder(json.JSONDecoder):
        def __init__(self, *args, **kwargs):
            super(CustomJSONDecoder, self).__init__(
                object_hook=self.object_hook, *args, **kwargs
            )

        def object_hook(self, dct):
            for coding_method in coding_methods:
                result = coding_method.object_hook(dct)
                if result is not dct:
                    return result
            return dct

    return (CustomJSONEncoder, CustomJSONDecoder)


class JSONCodec(object):
    """Custom JSON encoding and decoding for a specific class

    Parameters
    ----------
    cls : class
        Class for this codec. Assumes that all subclasses should be treated
        the same way. Can be `None` if `is_my_obj` and `is_my_class` are
        given.
    to_dict : callable
        method that converts the object to a dictionary
    from_dict : callable
        method that restores the object based on the dictionary made by
        to_dict
    is_my_obj : callable
        (Optional) Method to determine whether the input object should be
        treated by this encoder. Default behavior is to use
        ``isinstance(cls)``, and to create a dict that also includes the
        class name and the module of the object.
    is_my_dict : callable
        (Optional) Method to determine whether the input dictionary should
        be treated by this decoder. Default behavior assumes usage of the
        default ``is_my_obj``.
    """
    def __init__(self, cls, to_dict, from_dict, is_my_obj=None,
                 is_my_dict=None):
        self.cls = cls
        self.to_dict = to_dict
        self.from_dict = from_dict
        self.is_my_obj = none_to_default(is_my_obj, self._is_my_obj)
        self.is_my_dict = none_to_default(is_my_dict, self._is_my_dict)

    def _is_my_dict(self, dct):
        is_custom = '__class__' in dct and '__module__' in dct
        if is_custom:
            return (dct['__class__'] == self.cls.__name__
                    and dct['__module__'] == self.cls.__module__)

    def _is_my_obj(self, obj):
        return isinstance(obj, self.cls)

    def default(self, obj):
        if self.is_my_obj(obj):
            dct = {}
            if self.cls:
                dct.update({'__class__': self.cls.__name__,
                            '__module__': self.cls.__module__})
            # we let the object override __class__ and __module__ if needed
            dct.update(self.to_dict(obj))
            return dct
        return obj

    def object_hook(self, dct):
        if self.is_my_dict(dct):
            obj = self.from_dict(dct)
            return obj
        return dct

def iterable_to_dict(obj):
    return {'as_list': list(obj)}

def iterable_from_dict(dct, iterable_type):
    return iterable_type(dct['as_list'])

tuple_codec = JSONCodec(tuple, iterable_to_dict,
                        functools.partial(iterable_from_dict,
                                          iterable_type=tuple))
set_codec = JSONCodec(set, iterable_to_dict,
                        functools.partial(iterable_from_dict,
                                          iterable_type=set))

def bytes_to_dict(obj):
    return {'bytes': obj.decode('latin-1')}

def bytes_from_dict(dct):
    return dct['bytes'].encode('latin-1')

bytes_codec = JSONCodec(bytes, bytes_to_dict, bytes_from_dict)


import numpy as np
def numpy_to_dict(obj):
    return {'shape': obj.shape,
            'dtype': str(obj.dtype),
            'string': obj.tobytes()}

def numpy_from_dict(dct):
    arr = np.frombuffer(dct['string'], dtype=np.dtype(dct['dtype']))
    return arr.reshape(dct['shape'])

numpy_codec = JSONCodec(np.ndarray, numpy_to_dict, numpy_from_dict)


def uuid_object_to_dict(obj):
    dct = obj.to_dict()
    dct = replace_uuid(dct, uuid_encoding=encode_uuid)
    dct.update({'__class__': obj.__class__.__name__,
                '__module__': obj.__class__.__module__})
    name = getattr(obj, '_name', None)
    if name and 'name' not in dct:
        dct['name'] = name
    return dct

# we ignore all dicts on reserialization because we need to use the custom
# restoration process for that; so is_my_dict returns False
uuid_object_codec = JSONCodec(cls=None,
                              to_dict=uuid_object_to_dict,
                              from_dict=lambda x: x,  # never called
                              is_my_obj=has_uuid,
                              is_my_dict=lambda x: False)

DEFAULT_CODECS = [
    uuid_object_codec,
    bytes_codec,
    tuple_codec,
    set_codec,
    numpy_codec,
]
