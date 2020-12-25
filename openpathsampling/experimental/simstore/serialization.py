import numpy as np
from .my_types import parse_ndarray_type
from . import serialization_helpers as serialization
from . import attribute_handlers
import json

import logging
logger = logging.getLogger(__name__)


def load_list_uuid(json_str, cache_list):
    uuid_list = json.loads(json_str)
    if uuid_list is None:
        return uuid_list
    uuid_list = [serialization.decode_uuid(u) for u in uuid_list]
    return [serialization.search_caches(uuid, cache_list)
            for uuid in uuid_list]



class SchemaDeserializer(object):
    default_handlers = {
        'lazy': serialization.search_caches,
        'uuid': serialization.search_caches,
        'list_uuid': load_list_uuid,
    }

    def __init__(self, schema, table, cls, handlers):
        self.schema = schema
        self.table = table
        if table is not None:
            self.entries = schema[table]
        else:
            self.entries = []
        self.cls = cls
        self.handler_factories = handlers
        self.attribute_handlers = self.init_attribute_handlers()

    # TODO: move this external
    # @staticmethod
    # def make_numpy_handler(dtype, shape):
        # return lambda data, _: np.fromstring(data, dtype=dtype).reshape(shape)

    def get_handler_from_factories(self, type_name):
        for factory in self.handler_factories:
            handler = factory.from_type_string(type_name)
            if handler is not None:
                return handler.deserialize

    def init_attribute_handlers(self):
        attribute_handlers = {}
        for (attr, type_name) in self.entries:
            handler = None
            if type_name in self.default_handlers:
                handler = self.default_handlers[type_name]
            else:
                handler = self.get_handler_from_factories(type_name)
                # as_ndarray = parse_ndarray_type(type_name)
                # if as_ndarray:
                    # (dtype, shape) = as_ndarray
                    # handler = self.make_numpy_handler(dtype, shape)
            if handler:
                attribute_handlers[attr] = handler
        return attribute_handlers

    def make_dct(self, table_dct, cache_list):
        for attr in self.attribute_handlers:
            table_dct[attr] = self.attribute_handlers[attr](table_dct[attr],
                                                            cache_list)
        return table_dct

    def __call__(self, uuid, table_dct, cache_list):
        dct = self.make_dct(table_dct, cache_list)
        # if 'uuid' in dct:
            # del dct['uuid']
        obj = self.cls.from_dict(dct)
        serialization.set_uuid(obj, uuid)
        return obj


class ToDictSerializer(SchemaDeserializer):
    default_handlers = {
        'uuid': serialization.get_uuid,
        'lazy': serialization.get_uuid,
        'json': serialization.to_bare_json,
        'json_obj': serialization.to_json_obj,
        'list_uuid': serialization.to_bare_json
    }

    # TODO: I think this class can be basically remobed; need to transfer
    # inherited func and attribs to SchemaSerializer

    def get_handler_from_factories(self, type_name):
        for factory in self.handler_factories:
            handler = factory.from_type_string(type_name)
            if handler is not None:
                return handler.serialize

    def __call__(self, obj):
        dct = obj.to_dict()
        for attr in self.attribute_handlers:
            dct[attr] = self.attribute_handlers[attr](dct[attr])
        dct = serialization.replace_uuid(dct,
                                         uuid_encoding=lambda x: x)
        dct.update({'uuid': serialization.get_uuid(obj)})
        return dct


class SchemaSerializer(ToDictSerializer):
    def __call__(self, obj):
        dct = {attr: getattr(obj, attr)
               for (attr, type_name) in self.entries}
        replace = {attr: handler(dct[attr])
                   for (attr, handler) in self.attribute_handlers.items()}
        replace = serialization.replace_uuid(replace,
                                             uuid_encoding=lambda x: x)
        dct.update(replace)
        dct.update({'uuid': serialization.get_uuid(obj)})
        return dct
