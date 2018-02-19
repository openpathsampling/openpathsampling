from serialization_helpers import *
import numpy as np
import pytest

class MockUUIDObject(object):
    def __init__(self, uuid, normal_attr=None, obj_attr=None,
                 list_attr=None, dict_attr=None):
        self.__uuid__ = uuid
        self.dict_attr = dict_attr
        self.list_attr = list_attr
        self.obj_attr = obj_attr
        self.normal_attr = normal_attr

    def to_dict(self):
        return {
            'obj_attr': self.obj_attr,
            'list_attr': self.list_attr,
            'dict_attr': self.dict_attr,
            'normal_attr': self.normal_attr
        }

    @classmethod
    def from_dict(cls, dct):
        # set UUID after
        return cls(uuid=None, **dct)

def create_test_objects():
    obj_int = MockUUIDObject(uuid='int', normal_attr=5)
    obj_str = MockUUIDObject(uuid='str', normal_attr='foo')
    obj_np = MockUUIDObject(uuid='np', normal_attr=np.array([1.0, 2.0]))
    obj_obj = MockUUIDObject(uuid='obj', obj_attr=obj_int)
    obj_lst = MockUUIDObject(uuid='lst', list_attr=[obj_int, obj_str])
    obj_dct = MockUUIDObject(uuid='dct', dict_attr={'foo': obj_str,
                                                    obj_int: obj_np})
    obj_nest = MockUUIDObject(
        uuid='nest',
        dict_attr={'bar': [obj_str, {obj_int: [obj_np, obj_obj]}]}
    )
    obj_repeat = MockUUIDObject('rep', list_attr=[obj_int, [obj_int]])
    all_objects = {
        obj.__uuid__: obj
        for obj in [obj_int, obj_str, obj_np, obj_obj, obj_lst, obj_dct,
                    obj_nest, obj_repeat]
    }
    return all_objects

all_objects = create_test_objects()


@pytest.mark.parametrize('obj', list(all_objects.values()))
def test_has_uuid(obj):
    assert has_uuid(obj)

def test_has_uuid_no_uuid():
    assert not has_uuid(10)
    assert not has_uuid('foo')

@pytest.mark.parametrize('uuid,obj', list(all_objects.items()))
def test_get_uuid(uuid, obj):
    assert get_uuid(obj) == uuid

def test_get_uuid_none():
    assert get_uuid(None) == None

@pytest.mark.parametrize('obj,included_objs', [
    (all_objects['int'], []),
    (all_objects['str'], []),
    (all_objects['np'], []),
    (all_objects['obj'], [all_objects['int']]),
    (all_objects['lst'], [all_objects['int'], all_objects['str']]),
    (all_objects['dct'], [all_objects['str'], all_objects['int'],
                          all_objects['np']]),
    (all_objects['nest'], [all_objects['str'], all_objects['int'],
                           all_objects['np'], all_objects['obj']]),
    (all_objects['rep'], [all_objects['int']])
])
def test_get_all_uuids(obj, included_objs):
    expected = {o.__uuid__: o for o in included_objs}
    expected.update({obj.__uuid__: obj})
    assert get_all_uuids(obj) == expected

@pytest.mark.parametrize('obj,included_objs', [
    (all_objects['int'], []),
    (all_objects['str'], []),
    (all_objects['np'], []),
    (all_objects['obj'], []),
    (all_objects['lst'], [all_objects['str']]),
    (all_objects['dct'], [all_objects['str'], all_objects['np']]),
    (all_objects['nest'], [all_objects['str'], all_objects['np'],
                           all_objects['obj']]),
    (all_objects['rep'], [])
])
def test_get_all_uuids_with_known(obj, included_objs):
    expected = {o.__uuid__: o for o in included_objs}
    known_obj = all_objects['int']
    known_uuids = {known_obj.__uuid__: known_obj}
    if obj is not all_objects['int']:
        expected.update({obj.__uuid__: obj})
    assert get_all_uuids(obj, known_uuids=known_uuids) == expected

# def test_replace_uuid(obj, replace_dct):
    # pytest.skip()

def test_search_caches():
    pytest.skip()

def test_search_caches_missing():
    pytest.skip()

def test_seach_caches_missing_error():
    pytest.skip()

def test_from_dict_with_uuids():
    pytest.skip()


