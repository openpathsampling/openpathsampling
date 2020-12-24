from collections import namedtuple
import json
from .serialization_helpers import *
from .serialization_helpers import _uuids_from_table_row
import numpy as np
import pytest

from .test_utils import (toy_uuid_maker, uuid_encode, MockUUIDObject,
                         MockBackend, create_test_objects)

all_objects = create_test_objects()

def get_obj_uuid(name):
    return get_uuid(all_objects[name])


@pytest.mark.parametrize('obj', list(all_objects.values()))
def test_has_uuid(obj):
    assert has_uuid(obj)

def test_has_uuid_no_uuid():
    assert not has_uuid(10)
    assert not has_uuid('foo')

@pytest.mark.parametrize('name,obj', list(all_objects.items()))
def test_get_uuid(name, obj):
    assert get_uuid(obj) == str(int(hash(name)))

def test_get_uuid_none():
    assert get_uuid(None) is None

def test_get_uuid_error():
    with pytest.raises(AttributeError):
        get_uuid(10)

def test_set_uuid():
    obj = MockUUIDObject(name="test", normal_attr=10)
    fake_uuid = 100
    assert get_uuid(obj) != str(fake_uuid)
    set_uuid(obj, fake_uuid)
    assert get_uuid(obj) == str(fake_uuid)

def test_encode_uuid():
    obj = MockUUIDObject(name="test", normal_attr=10)
    set_uuid(obj, 100)
    assert encode_uuid("UUID(100)")

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
    expected = {str(o.__uuid__): o for o in included_objs}
    expected.update({str(obj.__uuid__): obj})
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
    expected = {str(o.__uuid__): o for o in included_objs}
    known_obj = all_objects['int']
    known_uuids = {str(known_obj.__uuid__): known_obj}
    if obj is not all_objects['int']:
        expected.update({str(obj.__uuid__): obj})
    assert get_all_uuids(obj, known_uuids=known_uuids) == expected

@pytest.mark.parametrize('obj,replace_dct', [
    (all_objects['int'], {'name': 'int', 'normal_attr': 5}),
    (all_objects['str'], {'name': 'str', 'normal_attr': 'foo'}),
    (all_objects['obj'], {'name': 'obj',
                          'obj_attr': uuid_encode('int')}),
    (all_objects['lst'], {'name': 'lst',
                          'list_attr': [uuid_encode('int'),
                                        uuid_encode('str')]}),
    (all_objects['dct'], {'name': 'dct',
                          'dict_attr': {
                              'foo': uuid_encode('str'),
                              uuid_encode('int'): uuid_encode('np')
                          }}),
    (all_objects['nest'], {
        'name': 'nest',
        'dict_attr': {
            'bar': [uuid_encode('str'),
                    {uuid_encode('int'): [uuid_encode('np'),
                                          uuid_encode('obj')]}]
        }}),
    (all_objects['rep'], {'name': 'rep',
                          'list_attr': [uuid_encode('int'),
                                        [uuid_encode('int')]]})
])
def test_replace_uuid(obj, replace_dct):
    after_replacement = {key: None for key in MockUUIDObject.attr_list}
    after_replacement.update(replace_dct)
    encoding = lambda x: "UUID(" + str(x) + ")"
    assert replace_uuid(obj.to_dict(), encoding) == after_replacement

def test_replace_uuid_ndarray():
    # can't use assert == with arrays
    after_replacement = {'name': 'np', 'dict_attr': None, 'lazy_attr': None,
                         'list_attr': None, 'obj_attr': None,
                         'normal_attr': np.array([1.0, 2.0])}
    encoding = lambda x: "UUID(" + str(x) + ")"
    result = replace_uuid(all_objects['np'].to_dict(), encoding)
    assert set(result.keys()) == set(after_replacement.keys())
    for key in after_replacement:
        if key == 'normal_attr':
            assert np.allclose(result[key], after_replacement[key])
        else:
            assert result[key] == after_replacement[key]

@pytest.fixture
def cache_list():
    make_cache = lambda keys: {get_uuid(all_objects[key]): all_objects[key]
                               for key in keys}
    cache_1 = make_cache(['int', 'str'])
    cache_2 = make_cache(['obj'])
    return [cache_1, cache_2]

def test_search_caches(cache_list):
    for key in ['int', 'str', 'obj']:
        uuid = get_uuid(all_objects[key])
        assert search_caches(uuid, cache_list) == all_objects[key]

    # seearch in a single cache (listify)
    for key in ['int', 'str']:
        uuid = get_uuid(all_objects[key])
        assert search_caches(uuid, cache_list[0]) == all_objects[key]

def test_search_caches_missing(cache_list):
    assert search_caches("foo", cache_list, raise_error=False) is None
    uuid = get_uuid(all_objects['obj'])
    assert search_caches(uuid, cache_list[0], raise_error=False) is None

def test_seach_caches_missing_error(cache_list):
    with pytest.raises(KeyError):
        search_caches("foo", cache_list)

    with pytest.raises(KeyError):
        uuid = get_uuid(all_objects['obj'])
        search_caches(uuid, cache_list[0])

def test_from_dict_with_uuids(cache_list):
    # this one only uses lst
    # this matches test_replace_uuid
    lst_dict = {key: None for key in MockUUIDObject.attr_list}
    lst_dict.update({'name': 'lst', 'list_attr': [uuid_encode('int'),
                                                  uuid_encode('str')]})

    lst_obj = from_dict_with_uuids(lst_dict, cache_list)
    expected = all_objects['lst']
    assert lst_obj == expected.to_dict()

def test_uuids_from_table_row():
    TableRow = namedtuple("TableRow",
                          ['uuid', 'idx'] + MockUUIDObject.attr_list)
    row = TableRow(name="row",
                   dict_attr=None,
                   list_attr=json.dumps([uuid_encode('int'),
                                         uuid_encode('str')]),
                   obj_attr=str(toy_uuid_maker('obj')),
                   lazy_attr=str(toy_uuid_maker('nest')),
                   normal_attr="bar",
                   uuid=str(toy_uuid_maker('row')),
                   idx=1)
    entries = MockUUIDObject.schema
    uuids, lazy, deps = _uuids_from_table_row(row, entries)

    assert lazy == {str(toy_uuid_maker('nest'))}
    # TODO: do we want to allow None in the UUID list? Comes from dict_attr
    assert set(uuids) == {str(toy_uuid_maker('int')),
                          str(toy_uuid_maker('str')),
                          str(toy_uuid_maker('obj')), None}
    assert len(deps) == 1
    assert set(deps[str(toy_uuid_maker('row'))]) == \
            {str(toy_uuid_maker('int')), str(toy_uuid_maker('str')),
             str(toy_uuid_maker('obj')), str(toy_uuid_maker('nest'))}

def test_schema_find_uuids():
    test_objs = create_test_objects()
    test_objs['lazy'] = test_objs['obj']
    schemas = {
        'int': [('normal_attr', 'int')],
        'str': [('normal_attr', 'str')],
        'lst': [('list_attr', 'list_uuid')],
        'obj': [('obj_attr', 'uuid')],
        'lazy': [('obj_attr', 'lazy')]
    }
    expected_newobjs = {
        'int': [],
        'str': [],
        'lst': [test_objs['int'], test_objs['str']],
        'obj': [test_objs['int']],
        'lazy': [test_objs['int']]
    }
    for (key, schema_entries) in schemas.items():
        obj = test_objs[key]
        schema_find_uuids = SchemaFindUUIDs(schema_entries)
        uuids, new_objs = schema_find_uuids(test_objs[key], cache_list=[])
        assert uuids == {get_uuid(obj): obj}
        assert new_objs == expected_newobjs[key]

def test_get_all_uuids_loading():
    backend = MockBackend()
    schema = {
        'sims': [('json', 'json_obj'), ('class_idx', 'int')],
        'ints': [('normal_attr', 'int')],
        'objs': [('obj_attr', 'uuid')],
    }
    uuid_list = [get_obj_uuid('dct2'), get_obj_uuid('obj')]
    all_table_rows, lazy, dependencies, uuid_to_table = \
            get_all_uuids_loading(uuid_list, backend, schema)

    expected_dependencies = {
        get_obj_uuid('dct2'): {get_obj_uuid('str'), get_obj_uuid('int')},
        get_obj_uuid('obj'): {get_obj_uuid('int')},
        get_obj_uuid('str'): set([]),
        get_obj_uuid('int'): set([])
    }
    expected_uuid_to_table = {
        get_obj_uuid('dct2'): 'sims',
        get_obj_uuid('str'): 'sims',
        get_obj_uuid('int'): 'ints',
        get_obj_uuid('obj'): 'objs'
    }
    assert len(all_table_rows) == 4
    uuid_to_table_row = {row.uuid: row for row in all_table_rows}
    # we force the same index as the one with the same UUID, since we can't
    # know the index in advance
    expected_table_rows = []
    for obj in ['dct2', 'str', 'int', 'obj']:
        uuid = get_obj_uuid(obj)
        row_type = backend.row_types[expected_uuid_to_table[uuid]]
        # if this raises a KeyError, then the UUID we expected wasn't found
        idx = uuid_to_table_row[uuid].idx
        table_data = backend.table_data[obj]
        row_dict = {k: v for k, v in table_data.items()
                    if k not in ['obj', 'table_name']}
        row_dict.update({'idx': idx, 'uuid': uuid})
        row = row_type(**row_dict)
        expected_table_rows.append(row)

    assert set(all_table_rows) == set(expected_table_rows)
    assert lazy == set([])
    assert dependencies == expected_dependencies
    assert uuid_to_table == expected_uuid_to_table

@pytest.mark.parametrize('lazy_allowed', [True, False])
def test_get_all_uuids_loading_lazy_allowed(lazy_allowed):
    ...
    pytest.skip()

def test_get_all_uuids_loading_with_existing():
    pytest.skip()

def test_dependency_dag():
    # check that we have the expected nodes, edges
    pytest.skip()

def test_dependency_dag_with_initial_dag():
    pytest.skip()

def test_get_reload_order():
    # check order, including no-dep orders
    pytest.skip()
