import numpy as np
from my_types import parse_ndarray_type, ndarray_re
import serialization_helpers as serialization
from tools import is_mappable

def make_lazy_class(cls_):
    # this is to mix-in inheritence
    class LazyClass(LazyLoader, cls_):
        pass
    return LazyClass

class LazyLoader(object):
    def __init__(self, uuid, class_, storage):
        self.__uuid__ = uuid
        self.storage = storage
        self.class_= class_
        self._loaded_object = None

    def load(self):
        if self._loaded_object is None:
            self._loaded_object = self.storage.load(self.__uuid__, lazy=False)
        if self._loaded_object is None:
            raise RuntimeError("UUID not found in storage: " +
                               str(self.__uuid__))
        return self._loaded_object

    def __getattr__(self, attr):
        return getattr(self.load(), attr)

    def __iter__(self):
        loaded = self.load()
        if not hasattr(loaded, '__iter__'):
            raise TypeError() # TODO: message
        # TODO: load all objects in the list
        return loaded.__iter__

    # TODO: how to handle isinstance? each lazy-loading class needs a sub
    def repr(self):
        return ("<LazyLoader for " + str(self.class_.__name__) + " UUID "
                + str(self.__uuid__) + ">")

class Serialization(object):
    builtin_types = ['int', 'float', 'str']
    uuid_types = ['uuid', 'list_uuid', 'lazy']
    # TODO: to_json here might not quite be correct; need to_bare_json?
    def __init__(self, storage):
        self.storage = storage
        self.cache = self.storage.cache
        self.attribute_serializers = {
            'uuid': serialization.get_uuid,
            'lazy': serialization.get_uuid,
            'json': serialization.to_bare_json,
            'list_uuid': serialization.to_bare_json
        }

        self.attribute_deserializers = {
            'uuid': serialization.from_json_obj,
            'lazy': self.make_lazy,
            'json': serialization.from_bare_json,
            'list_uuid': serialization.from_bare_json
        }
        self.schema = {}
        self.table_to_class = {}
        self._ser_dict = {}
        self._deser_dict = {}
        self._lazy_classes = {}

    def make_lazy(self, cls, uuid):
        if cls not in self._lazy_classes:
            self._lazy_classes[cls] = make_lazy_class(cls)
        return self._lazy_classes[cls](uuid=uuid,
                                       class_=cls,
                                       storage=self.storage)

    def register_serialization(self, schema, class_info):
        for table in schema:
            if class_info[table].serializer:
                self._ser_dict[table] = class_info[table].serializer
            else:
                self._ser_dict[table] = \
                        self.default_serializer_dict(schema[table])

            if class_info[table].deserializer:
                self._deser_dict[table] = class_info[table].deserializer
            else:
                self._deser_dict[table] = \
                        self.default_deserializer_dict(schema[table])

            self.table_to_class.update({table: class_info[table].cls})
            self.schema.update(schema)

    def attribute_serializer(self, type_name):
        if type_name in self.attribute_serializers:
            return self.attribute_serializers[type_name]
        if ndarray_re.match(type_name):
            return lambda arr: arr.tostring()
        else:
            raise TypeError("Unknown type for serialization: " + type_name)

    def attribute_deserializer(self, type_name):
        if type_name in self.attribute_deserializers:
            return self.attribute_deserializers[type_name]
        as_ndarray = parse_ndarray_type(type_name)
        if as_ndarray:
            (dtype, shape) = as_ndarray
            return lambda data: \
                    np.fromstring(data, dtype=dtype).reshape(shape)

    def _serialization_dict(self, attribute_handler, table_description):
        dct = {}
        for (attr, type_name) in table_description:
            if type_name in self.builtin_types:
                dct[attr] = None
            else:
                dct[attr] = attribute_handler(type_name)
        return dct

    def default_serializer_dict(self, table_description):
        return self._serialization_dict(self.attribute_serializer,
                                        table_description)

    def default_deserializer_dict(self, table_description):
        return self._serialization_dict(self.attribute_deserializer,
                                        table_description)

    @property
    def serialize(self):
        default_tables = set(tab for (tab, func) in self._ser_dict.items()
                             if is_mappable(func))
        non_default_tables = set(self._ser_dict.keys()) - default_tables
        results = {table: self._ser_dict[table]
                   for table in non_default_tables}
        results.update({
            table:
            lambda obj, table=table: self.default_serialize(obj, table)
            for table in default_tables
        })
        return results

    def default_serialize(self, obj, table):
        all_objects = []
        dct = {'uuid': serialization.get_uuid(obj)}
        serializer_dict = self._ser_dict[table]
        for (attr, type_name) in self.schema[table]:
            attr_obj = getattr(obj, attr)
            if attr not in serializer_dict:
                dct[attr] = attr_obj  # built-in types
            else:
                dct[attr] = serializer_dict[attr](attr_obj)
        return dct

    def default_deserialize(self, uuid, dct, table):
        # this must be called in the DAG order; otherwise error
        deserializer_dict = self._deser_dict[table]
        cls = self.table_to_class[table]
        for (attr, type_name) in self.schema[table]:
            if type_name == 'lazy':
                dct[attr] = self.make_lazy(cls, uuid)
            elif type_name in ['uuid', 'list_uuid']:
                dct[attr] = deserializer_dict[attr](dct[attr], self.cache)
            else:
                dct[attr] = deserializer_dict[attr](dct[attr])
        obj = cls.from_dict(dct)
        obj.__uuid__ = uuid
        return obj

