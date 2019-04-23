"""Utilities for tests: mocks, etc.
"""
from collections import namedtuple
import pytest
import numpy as np

from .serialization_helpers import get_uuid


def toy_uuid_maker(name):
    return int(hash(name))

def uuid_encode(name):
    return "UUID(" + str(toy_uuid_maker(name)) + ")"


class AbstractMockUUIDObject(object):
    def __init__(self, *args, **kwargs):
        keywords = dict(zip(self.attr_list, args))
        check = [attr in self.attr_list and attr not in keywords
                 for attr in kwargs]
        if all(check):
            keywords.update(kwargs)
        else:
            raise Exception("Something bad happened in setup")

        # all defaults are None
        keywords.update({attr: None for attr in self.attr_list
                         if attr not in keywords})

        self.__uuid__ = toy_uuid_maker(keywords['name'])

        for attr, value in keywords.items():
            setattr(self, attr, value)

    def to_dict(self):
        return {attr: getattr(self, attr) for attr in self.attr_list}

    @classmethod
    def from_dict(cls, dct):
        # set UUID after
        return cls(name=None, **dct)

class MockUUIDObject(AbstractMockUUIDObject):
    attr_list = ['name', 'normal_attr', 'obj_attr', 'list_attr',
                 'dict_attr', 'lazy_attr']
    schema = [('dict_attr', 'uuid'), ('list_attr', 'list_uuid'),
              ('obj_attr', 'uuid'), ('lazy_attr', 'lazy'),
              ('normal_attr', 'str')]

class MockSimulationObject(AbstractMockUUIDObject):
    attr_list = ['name', 'normal_attr']
    # no schema; use this for simulation objects

class ExtraMockDataObject(AbstractMockUUIDObject):
    attr_list = ['name', 'str_attr']
    schema = [('str_attr', 'str')]


class MockBackend(object):
    def _table_data_for_object(self, obj, table_name, **kwargs):
        schema_entries = self.schema[table_name]
        row_type = self.row_types[table_name]
        uuid = get_uuid(obj)
        table_idx = self.table_names.index(table_name)
        idx = len(self.tables[table_idx])
        row_kwargs = {'uuid': uuid, 'idx': idx}
        row_kwargs.update(kwargs)
        row = row_type(**row_kwargs)
        uuid_row = self.row_types['uuids'](uuid=uuid, table=table_idx,
                                           idx=idx)
        return uuid_row, row

    def __init__(self):
        self.schema = {
            'objs': [('obj_attr', 'uuid')],
            'ints': [('normal_attr', 'int')],
            'sims': [('json', 'json_obj'), ('class_idx', 'int')]
        }
        self.row_types = {
            'uuids': namedtuple('UUIDsRow', ['uuid', 'table', 'idx']),
            'objs': namedtuple('ObjRow', ['uuid', 'idx', 'obj_attr']),
            'ints': namedtuple('IntRow', ['uuid', 'idx', 'normal_attr']),
            'sims': namedtuple('SimRow',
                               ['uuid', 'idx', 'json', 'class_idx']),
            'clss': namedtuple('ClsRow', ['idx', 'module', 'cls'])
        }
        self.obj_to_table = {
            'int': 'ints',
            'obj': 'objs',
            'dct': 'sims',
            'str': 'sims'
        }
        self.table_names = ['uuids', 'clss', 'sims', 'objs', 'ints']

        self.uuid_table = {}
        self.tables = [[] for table in self.table_names]
        uuid_row, table_row = self._table_data_for_object(
            obj=all_objects['int'],
            table_name='ints',
            normal_attr=5
        )
        self.uuid_table[uuid_row.uuid] = uuid_row
        self.tables[uuid_row.table].append(table_row)

    def load_uuids_table(self, new_uuids):
        return [self.uuid_table[uuid] for uuid in new_uuids]

    def load_table_data(self, uuid_rows):
        return [self.tables[row.table][row.idx] for row in uuid_rows]

    def uuid_row_to_table_name(self, row):
        return self.table_names[row.table]

def create_test_objects():
    obj_int = MockUUIDObject(name='int', normal_attr=5)
    obj_str = MockUUIDObject(name='str', normal_attr='foo')
    obj_np = MockUUIDObject(name='np', normal_attr=np.array([1.0, 2.0]))
    obj_obj = MockUUIDObject(name='obj', obj_attr=obj_int)
    obj_lst = MockUUIDObject(name='lst', list_attr=[obj_int, obj_str])
    obj_dct = MockUUIDObject(name='dct', dict_attr={'foo': obj_str,
                                                    obj_int: obj_np})
    obj_nest = MockUUIDObject(
        name='nest',
        dict_attr={'bar': [obj_str, {obj_int: [obj_np, obj_obj]}]}
    )
    obj_repeat = MockUUIDObject('rep', list_attr=[obj_int, [obj_int]])
    all_objects = {
        obj.name : obj
        for obj in [obj_int, obj_str, obj_np, obj_obj, obj_lst, obj_dct,
                    obj_nest, obj_repeat]
    }
    return all_objects


all_objects = create_test_objects()


class TestMockBackend(object):
    def setup(self):
        self.backend = MockBackend()

    @pytest.mark.parametrize(('obj_name', 'table_idx', 'idx'),
                             [('int', 4, 0)])
    def test_load_uuids_table(self, obj_name, table_idx, idx):
        # TODO: add more examples to test
        uuid = get_uuid(all_objects[obj_name])
        assert self.backend.load_uuids_table([uuid]) == \
                [(uuid, table_idx, idx)]

    def test_load_table_data(self):
        pass

    def test_uuid_row_to_table_name(self):
        pass
