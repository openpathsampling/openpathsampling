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
        self.info_default = ClassInfo(
            'simulation_objects', cls=MockSimulationObject,
            # TODO: serializer=ser,
            # TODO: deserializer=deser,
            find_uuids=default_find_uuids
        )
        self.info_mock = ClassInfo(table='mock', cls=MockUUIDObject)
        schema = {'mock': MockUUIDObject.schema}
        self.serialization_schema = SerializationSchema(
            default_info=self.info_default,
            schema=schema,
            class_info_list=[self.info_mock]
        )

    def test_init(self):
        serialization = self.serialization_schema

        schema = {'mock': MockUUIDObject.schema}
        assert serialization.schema == schema

        class_info_set = set([self.info_mock, self.info_default])
        assert set(serialization.class_info_list) == class_info_set

        table_to_info = {'mock': self.info_mock,
                         'simulation_objects': self.info_default}
        assert serialization.table_to_info == table_to_info

        lookup_to_info = {MockUUIDObject: self.info_mock,
                          MockSimulationObject: self.info_default}
        assert serialization.lookup_to_info == lookup_to_info

    def test_tables(self):
        expected = set(['simulation_objects', 'mock'])
        assert set(self.serialization_schema.tables) == expected

    def test_register_info(self):
        serialization = self.serialization_schema
        new_class_info = ClassInfo(table='extra', cls=ExtraMockDataObject)

        assert new_class_info.serializer is None
        assert new_class_info.deserializer is None
        assert new_class_info.find_uuids is None
        new_schema = {'extra': ExtraMockDataObject.schema}
        serialization.register_info([new_class_info], schema=new_schema)

        # check that default behavior was set on the class_info
        assert new_class_info.serializer is not None
        assert new_class_info.deserializer is not None
        assert new_class_info.find_uuids is not None

        # check that the serialization schema was correctly updated
        schema = {'mock': MockUUIDObject.schema,
                  'extra': ExtraMockDataObject.schema}
        assert serialization.schema == schema

        class_info_set = set([self.info_mock, self.info_default,
                              new_class_info])
        assert set(serialization.class_info_list) == class_info_set

        table_to_info = {'mock': self.info_mock,
                         'simulation_objects': self.info_default,
                         'extra': new_class_info}
        assert serialization.table_to_info == table_to_info

        lookup_to_info = {MockUUIDObject: self.info_mock,
                          MockSimulationObject: self.info_default,
                          ExtraMockDataObject: new_class_info}
        assert serialization.lookup_to_info == lookup_to_info

    def test_lookup_key(self):
        serialization = self.serialization_schema
        assert serialization.lookup_key(self.data_obj) == MockUUIDObject
        assert serialization.lookup_key(self.sim_obj) == MockSimulationObject
        assert serialization.lookup_key(5) == int

    @pytest.mark.parametrize('instance', ['int', 'data', 'sim', 'unknown'])
    def test_info_from_instance(self, instance):
        serialization = self.serialization_schema
        input_val, expected = {
            'int': (5, None),
            'data': (self.data_obj, self.info_mock),
            'sim': (self.sim_obj, self.info_default),
            'unknown': (self.extra_obj, None)
        }[instance]
        assert serialization.info_from_instance(input_val) == expected

    @pytest.mark.parametrize('in_type', ['string', 'data_obj', 'sim_obj'])
    def test_get_item(self, in_type):
        inp_val, expected = {
            'string': ('mock', self.info_mock),
            'data_obj': (self.data_obj, self.info_mock),
            'sim_obj': (self.sim_obj, self.info_default)
        }[in_type]
        assert self.serialization_schema[inp_val] == expected

    def test_get_item_key_error(self):
        with pytest.raises(KeyError):
            self.serialization_schema[self.extra_obj]


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

