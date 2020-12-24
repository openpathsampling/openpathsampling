"""Utilities for tests: mocks, etc.
"""
from collections import namedtuple
import pytest
import numpy as np
import json
import random

from .serialization_helpers import get_uuid, encode_uuid


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

        if 'name' in keywords:
            self.__uuid__ = toy_uuid_maker(keywords['name'])
            self.is_named = True
        else:
            self.__uuid__ = random.randint(0, 2**32)
            self.is_named = False

        for attr, value in keywords.items():
            setattr(self, attr, value)


    def to_dict(self):
        return {attr: getattr(self, attr) for attr in self.attr_list}

    @classmethod
    def from_dict(cls, dct):
        # set UUID after
        try:
            name = dct.pop('name')
        except KeyError:
            name = None
        return cls(name=name, **dct)

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


class UnnamedUUID(AbstractMockUUIDObject):
    attr_list = ['normal_attr']
    schema = [('normal_attr', 'int')]


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

        # this is where we add objects to the table
        self.table_data = {
            'int': {'obj': all_objects['int'], 'table_name': 'ints',
                    'normal_attr': 5},
            'obj': {'obj': all_objects['obj'], 'table_name': 'objs',
                    'obj_attr': get_uuid(all_objects['int'])},
            'str': {'obj': all_objects['str'], 'table_name': 'sims',
                    'json': json.dumps('foo'), 'class_idx': 1},
            'dct2': {'obj': all_objects['dct2'], 'table_name': 'sims',
                     'json': json.dumps({'dct_attr': {
                         'foo': encode_uuid(toy_uuid_maker('str')),
                         encode_uuid(toy_uuid_maker('int')): 5
                     }}),
                     'class_idx': 0}
        }
        for datum in self.table_data.values():
            uuid_row, table_row = self._table_data_for_object(**datum)
            self.uuid_table[uuid_row.uuid] = uuid_row
            self.tables[uuid_row.table].append(table_row)

    def load_uuids_table(self, new_uuids):
        return [self.uuid_table[uuid] for uuid in new_uuids]

    def load_table_data(self, uuid_rows):
        return [self.tables[row.table][row.idx] for row in uuid_rows]

    def uuid_row_to_table_name(self, row):
        return self.table_names[row.table]


class LoadingStorageMock(object):
    def __init__(self, uuid_dict):
        self.uuid_dict = uuid_dict

    def load(self, uuid_list, force=False):
        return [self.uuid_dict[uuid] for uuid in uuid_list]

def create_test_objects():
    obj_int = MockUUIDObject(name='int', normal_attr=5)
    obj_str = MockUUIDObject(name='str', normal_attr='foo')
    obj_np = MockUUIDObject(name='np', normal_attr=np.array([1.0, 2.0]))
    obj_obj = MockUUIDObject(name='obj', obj_attr=obj_int)
    obj_lst = MockUUIDObject(name='lst', list_attr=[obj_int, obj_str])
    obj_dct = MockUUIDObject(name='dct', dict_attr={'foo': obj_str,
                                                    obj_int: obj_np})
    obj_dct2 = MockUUIDObject(name='dct2', dict_attr={'foo': obj_str,
                                                      obj_int: 5})
    obj_nest = MockUUIDObject(
        name='nest',
        dict_attr={'bar': [obj_str, {obj_int: [obj_np, obj_obj]}]}
    )
    obj_repeat = MockUUIDObject('rep', list_attr=[obj_int, [obj_int]])
    all_objects = {
        obj.name : obj
        for obj in [obj_int, obj_str, obj_np, obj_obj, obj_lst, obj_dct,
                    obj_dct2, obj_nest, obj_repeat]
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

    @pytest.mark.parametrize('table,uuid,table_idx,idx,content', [
        ('ints', 'int', 4, 0, [5])
    ])
    def test_load_table_data(self, table, uuid, table_idx, idx, content):
        row_types = self.backend.row_types
        uuid = str(toy_uuid_maker(uuid))
        uuid_row = row_types['uuids'](uuid=uuid, table=table_idx, idx=idx)
        expected = row_types[table](uuid, idx, *content)
        assert self.backend.load_table_data([uuid_row]) == [expected]

    @pytest.mark.parametrize('name,idx', [
        ('sims', 2), ('objs', 3), ('ints', 4)
    ])
    def test_uuid_row_to_table_name(self, name, idx):
        uuid_row = self.backend.row_types['uuids']('foo', idx, 0)
        assert self.backend.uuid_row_to_table_name(uuid_row) == name
