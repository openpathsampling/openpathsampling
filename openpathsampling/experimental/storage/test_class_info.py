from .class_info import *
import pytest

from .test_utils import (
    all_objects, MockUUIDObject, MockSimulationObject, ExtraMockDataObject
)

from .serialization_helpers import default_find_uuids

class TestClassInfo(object):
    def test_init_default_setup(self):
        class_info = ClassInfo(
            table="mock",
            cls=MockUUIDObject
        )
        schema = {'mock': MockUUIDObject.schema}
        class_info.set_defaults(schema)
        assert class_info.table == "mock"
        assert class_info.cls == MockUUIDObject
        assert isinstance(class_info.serializer, SchemaSerializer)
        assert isinstance(class_info.deserializer, SchemaDeserializer)
        assert isinstance(class_info.find_uuids, SchemaFindUUIDs)
        assert class_info.lookup_result == MockUUIDObject


class TestSerializationSchema(object):
    def setup(self):
        self.data_obj = all_objects['int']  # doesn't really matter which
        # technically, the next two should behave identically, except
        # name/UUID; but as different types they can be serailized
        # differently (data object serialization or simulation object
        # serialization)
        self.sim_obj = MockSimulationObject(name='sim_obj',
                                            normal_attr='foo')
        self.extra_obj = ExtraMockDataObject(name='data', str_attr='foo')
        class_info_list = [ClassInfo(table='mock', cls=MockUUIDObject)]
        schema = {'mock': MockUUIDObject.schema}
        self.serialization_schema = SerializationSchema(
            default_info=ClassInfo(
                'simulation_objects', cls=MockSimulationObject,
                # TODO: serializer=ser,
                # TODO: deserializer=deser,
                find_uuids=default_find_uuids
            ),
            schema=schema,
            class_info_list = class_info_list
        )


    def test_init(self):
        pytest.skip()

    def test_tables(self):
        pytest.skip()

    def test_add_class_info(self):
        pytest.skip()

    def test_register_info(self):
        pytest.skip()

    def test_lookup_key(self):
        pytest.skip()

    def test_info_from_instance(self):
        pytest.skip()


class TestSerializationSchemaSpecial(object):
    def test_special_lookup_key(self):
        pytest.skip()

    def test_is_special(self):
        pytest.skip()

    def test_get_special(self):
        pytest.skip()

    def test_lookup_key_for_special(self):
        pytest.skip()

    def test_info_from_instance_for_special(self):
        pytest.skip()

    def test_add_missing_table_from_instance(self):
        pytest.skip()

