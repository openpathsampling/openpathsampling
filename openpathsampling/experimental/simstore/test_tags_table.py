import pytest
from .storage import GeneralStorage
from .sql_backend import SQLStorageBackend
from .class_info import ClassInfo, SerializationSchema
from .custom_json import JSONSerializerDeserializer, uuid_object_codec
from .serialization_helpers import default_find_uuids

from .attribute_handlers import DEFAULT_HANDLERS

from openpathsampling.netcdfplus import StorableObject

from .tags_table import *

class IntHolder(StorableObject):
    def __init__(self, value):
        super().__init__()
        self.value = value

class TestTagsTable(object):
    def setup(self):
        json_ser = JSONSerializerDeserializer([uuid_object_codec])
        # TODO: add tags to serialization schema
        serialization_schema = SerializationSchema(
            default_info=ClassInfo(
                table='simulation_objects',
                cls=StorableObject,
                serializer=json_ser.simobj_serializer,
                deserializer=json_ser.simobj_deserializer,
                find_uuids=default_find_uuids
            ),
            sfr_info=None,
            schema={'int_holders': [('value', 'int')],
                    '_tagged_content': [('json', 'json_obj')]
                   },
            class_info_list=[
                ClassInfo(table='int_holders', cls=IntHolder),
                ClassInfo(table='_tagged_content', cls=TaggedObject,
                          serializer=json_ser.simobj_serializer,
                          deserializer=json_ser.simobj_deserializer,
                          find_uuids=default_find_uuids)
            ]
        )
        for info in serialization_schema.class_info_list:
            info.set_defaults(serialization_schema.schema, DEFAULT_HANDLERS)
        backend = SQLStorageBackend(":memory:", mode='w')
        self.storage = GeneralStorage(backend, serialization_schema,
                                      serialization_schema.schema)

        self.tags_table = TagsTable(self.storage)

    def test_init(self):
        assert self.tags_table.tagged_objects == {}

    def _set_info_in_storage(self):
        wrapped_holder = TaggedObject(IntHolder(10)).named('with-uuid')
        tagged_num = TaggedObject(5).named('no-uuid')
        self.storage.save([wrapped_holder, tagged_num])
        self.storage.backend.add_tag('tags', 'with-uuid',
                                     get_uuid(wrapped_holder))
        self.storage.backend.add_tag('tags', 'no-uuid',
                                     get_uuid(tagged_num))
        assert self.storage.backend.table_len('_tagged_content') == 2
        assert self.storage.backend.table_len('int_holders') == 1
        assert self.storage.backend.table_len('tags') == 2

    def test_load_tags(self):
        self._set_info_in_storage()
        tags_table = TagsTable(self.storage)
        assert len(tags_table.tagged_objects) == 2

    def test_len(self):
        assert len(self.tags_table) == 0
        self._set_info_in_storage()
        tags_table = TagsTable(self.storage)
        assert len(tags_table) == 2

    def test_getitem(self):
        self._set_info_in_storage()
        tags_table = TagsTable(self.storage)
        assert tags_table['no-uuid'] == 5
        with_uuid = tags_table['with-uuid']
        assert isinstance(with_uuid, IntHolder)
        assert with_uuid.value == 10

    def test_setitem(self):
        self.tags_table['no-uuid'] = 5  # JSON-serializable data
        self.tags_table['with-uuid'] = IntHolder(10)  # UUID objects
        assert len(self.tags_table.tagged_objects) == 2
        assert self.storage.backend.table_len('_tagged_content') == 2
        assert self.storage.backend.table_len('int_holders') == 1
        assert self.storage.backend.table_len('tags') == 2

    def test_setitem_name_error(self):
        self.tags_table['foo'] = 5
        with pytest.raises(RuntimeError):
            self.tags_table['foo'] = 6

    def test_set_get_cycle(self):
        wrapped = IntHolder(10)
        self.tags_table['no-uuid'] = 5
        self.tags_table['with-uuid'] = wrapped
        assert self.tags_table['no-uuid'] == 5
        assert self.tags_table['with-uuid'] == wrapped
