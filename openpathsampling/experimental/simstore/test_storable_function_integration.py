"""Integration tests for storable functions.

This includes more complicated tests of storable functions to ensure that
certain behaviors are preserved when interacting with storage.
"""

import pytest
import os

from openpathsampling.experimental.simstore import (
    StorableFunction, GeneralStorage, MemoryStorageBackend,
    SQLStorageBackend, ClassInfo, ClassInfoContainer,
    storable_function_find_uuids, StorableFunctionResults,
    CallableCodec
)
from openpathsampling.experimental.simstore.custom_json import (
    JSONSerializerDeserializer, DEFAULT_CODECS
)

from openpathsampling.experimental.simstore.serialization_helpers import \
        default_find_uuids

from openpathsampling.experimental.simstore.uuids import get_uuid

from openpathsampling.experimental.simstore.attribute_handlers import \
    DEFAULT_HANDLERS


from openpathsampling.netcdfplus import StorableNamedObject

class InputObj(StorableNamedObject):
    pass  # just a thing that has a UUID

_schema = {
    'input_objs': [],
    'storable_functions': [('json', 'json_obj')],
    'simulation_objects': [('json', 'json_obj')],
}

_json = JSONSerializerDeserializer(DEFAULT_CODECS + [CallableCodec()])

_serialization = ClassInfoContainer(
    default_info=ClassInfo(
        table='simulation_objects',
        cls=StorableNamedObject,
        serializer=_json.simobj_serializer,
        deserializer=_json.simobj_deserializer,
        find_uuids=default_find_uuids
    ),
    sfr_info=ClassInfo(
        table="function_results",
        cls=StorableFunctionResults,
        serializer=StorableFunctionResults.to_dict,
        deserializer=lambda x: x,  # not used
        find_uuids=default_find_uuids
    ),
    schema=_schema,
    class_info_list=[
        ClassInfo(table='input_objs', cls=InputObj),
        ClassInfo(table='storable_functions', cls=StorableFunction,
                  find_uuids=storable_function_find_uuids,
                  serializer=_json.simobj_serializer,
                  deserializer=_json.simobj_deserializer),
    ]
)

for info in _serialization.class_info_list:
    info.set_defaults(_schema, DEFAULT_HANDLERS)


class Mock(StorableNamedObject):
    # can't serialize mock.Mock, so we make our own
    def __init__(self, return_value):
        super().__init__()
        self.call_count = 0
        self.return_value = return_value

    def __call__(self, obj):
        self.call_count += 1
        return self.return_value


def test_setup():
    # test that we can store and load at all
    backend = MemoryStorageBackend()
    storage = GeneralStorage(backend=backend,
                             class_info=_serialization,
                             schema=_schema)
    obj = InputObj()
    assert storage.backend.table_len('input_objs') == 0
    storage.save(obj)
    assert storage.backend.table_len('input_objs') == 1
    func = Mock(return_value='foo')
    sf = StorableFunction(func).named('foo-return')
    sf_uuid = get_uuid(sf)
    assert not storage.backend.has_table(sf_uuid)
    assert len(sf.local_cache) == 0
    assert sf(obj) == 'foo'
    assert len(sf.local_cache) == 1
    storage.save(sf)
    assert storage.backend.has_table(sf_uuid)


def test_multiple_storage(tmpdir):
    # TODO: This should probably be refactored into several individual tests
    # that are each a little smaller. However, the main goal is to ensure
    # that this functionality gets tested
    inp1 = InputObj()
    inp2 = InputObj()
    storage = GeneralStorage(
        backend=SQLStorageBackend(tmpdir.join("st1.db"), mode='w'),
        class_info=_serialization,
        schema=_schema
    )
    func = Mock(return_value='foo')
    sf = StorableFunction(func).named('foo-return')
    out1 = sf(inp1)
    storage.save(sf)
    GeneralStorage._known_storages = {}
    storage = GeneralStorage(
        backend=SQLStorageBackend(tmpdir.join("st1.db"), mode='r'),
        class_info=_serialization,
        schema=_schema
    )
    sf = storage.storable_functions[0]
    assert sf.name == 'foo-return'
    assert len(sf.local_cache) == 0
    assert sf.func.call_count == 1

    # load result from storage: don't call func; do cache result
    _ = sf(inp1)
    assert sf.func.call_count == 1
    assert len(sf.local_cache) == 1

    # save to a new storage: this should save the cached result, too!
    _ = sf(inp2)
    assert sf.func.call_count == 2
    assert len(sf.local_cache) == 2
    st2 = GeneralStorage(
        backend=SQLStorageBackend(tmpdir.join("st2.db"), mode='w'),
        class_info=_serialization,
        schema=_schema
    )
    st2.save(sf)
    assert st2.backend.table_len(get_uuid(sf)) == 2

    # save one result to one storage, the other result to another storage;
    # can we get both results without calling the function?
    GeneralStorage._known_storages = {}
    storage = GeneralStorage(
        backend=SQLStorageBackend(tmpdir.join("st1.db"), mode='r'),
        class_info=_serialization,
        schema=_schema
    )
    sf = storage.storable_functions[0]
    assert sf.name == 'foo-return'
    assert len(sf.local_cache) == 0
    assert sf.func.call_count == 1
    st2 = GeneralStorage(
        backend=SQLStorageBackend(tmpdir.join("st2.db"), mode='w'),
        class_info=_serialization,
        schema=_schema
    )
    _ = sf(inp2)
    st2.save(sf)
    assert sf.func.call_count == 2
    res1 = sf(inp1)
    assert sf.func.call_count == 2
    res2 = sf(inp2)
    assert sf.func.call_count == 2


class Container(StorableNamedObject):
    def __init__(self, cv):
        super().__init__()
        self.cv = cv

    def __call__(self, inp):
        return self.cv(inp)[0]

@pytest.fixture
def inputs_and_func():
    inp1 = InputObj()
    inp2 = InputObj()
    func = Mock(return_value='foo')
    sf = StorableFunction(func).named('foo-return')
    container = Container(sf).named("container")
    return inp1, inp2, func, sf, container

def _store_via_container(backend, container, inp, call_count,
                         backend_count):
    storage = GeneralStorage(
        backend=backend,
        class_info=_serialization,
        schema=_schema
    )
    assert container(inp) == 'f'
    assert container.cv.func.call_count == call_count
    storage.save(container)
    sf_uuid = str(get_uuid(container.cv))
    assert storage.backend.has_table(sf_uuid)
    assert storage.backend.table_len(sf_uuid) == backend_count
    return storage

def _load_container(storage, sf):
    container_list = [obj for obj in storage.simulation_objects
                      if obj.name == 'container']
    assert len(container_list) == 1
    container = container_list[0]
    assert get_uuid(container.cv) == get_uuid(sf)
    assert container.cv is not sf  # should not be the same obj in memory
    assert container.cv.func.call_count == 1
    assert len(container.cv.local_cache) == 0
    return container

def test_storage_inside_other_object(tmpdir, inputs_and_func):
    inp1, inp2, func, sf, container = inputs_and_func
    be_1w = SQLStorageBackend(tmpdir.join("st1.db"), mode='w')
    _store_via_container(be_1w, container, inp1, call_count=1,
                         backend_count=1)


    GeneralStorage._known_storages = {}
    storage = GeneralStorage(
        backend=SQLStorageBackend(tmpdir.join("st1.db"), mode='a'),
        class_info=_serialization,
        schema=_schema
    )
    container = _load_container(storage, sf)

    # same store; store new result should work
    assert container(inp2) == 'f'
    assert container.cv.func.call_count == 2
    storage.save([container, inp2])
    sf_uuid = str(get_uuid(sf))
    assert storage.backend.table_len(sf_uuid) == 2

def test_multiple_storage_inside_other_object(tmpdir, inputs_and_func):
    inp1, inp2, func, sf, container = inputs_and_func
    sf_uuid = str(get_uuid(sf))
    be_1w = SQLStorageBackend(tmpdir.join("st1.db"), mode='w')
    _store_via_container(be_1w, container, inp1, call_count=1,
                         backend_count=1)

    # different store, first store one result then store the other
    GeneralStorage._known_storages = {}
    storage = GeneralStorage(
        backend=SQLStorageBackend(tmpdir.join("st1.db"), mode='r'),
        class_info=_serialization,
        schema=_schema
    )
    container = _load_container(storage, sf)

    st2 = GeneralStorage(
        backend=SQLStorageBackend(tmpdir.join('st2.db'), mode='w'),
        class_info=_serialization,
        schema=_schema
    )
    # store container before using it: should not add table
    st2.save(container)
    assert not st2.backend.has_table(sf_uuid)
    assert len(container.cv.local_cache) == 0
    assert container.cv.func.call_count == 1

    # now we use container on inp2; should call CV; should disk-cache
    assert container(inp2) == 'f'
    assert container.cv.func.call_count == 2
    assert len(container.cv.local_cache) == 1
    st2.save([container, inp2])
    assert st2.backend.has_table(sf_uuid)
    assert st2.backend.table_len(sf_uuid) == 1

    # now we calculate and save the value with inp1; should load from
    # existing storage (st1)
    assert container.cv.func.call_count == 2
    assert len(container.cv.local_cache) == 1
    assert container(inp1) == 'f'
    assert container.cv.func.call_count == 2
    assert len(container.cv.local_cache) == 2
    st2.save([container, inp1])
    assert st2.backend.table_len(sf_uuid) == 2

    # now clear the local cache and load inp2; should load from st2 without
    # calculating
    container.cv.local_cache.clear()
    assert container.cv.func.call_count == 2
    assert container(inp2) == 'f'
    assert container.cv.func.call_count == 2 

def test_closing(tmpdir, inputs_and_func):
    inp1, inp2, func, sf, container = inputs_and_func
    sf_uuid = str(get_uuid(sf))
    st1_fname = tmpdir.join("st1.db")
    be_1w = SQLStorageBackend(st1_fname, mode='w')
    storage = _store_via_container(be_1w, container, inp1, call_count=1,
                                   backend_count=1)

    assert sf.has_handler
    storage.close()
    assert not sf.has_handler

def test_storage_file_problem(tmpdir, inputs_and_func):
    inp1, inp2, func, sf, container = inputs_and_func
    sf_uuid = str(get_uuid(sf))
    st1_fname = tmpdir.join("st1.db")
    be_1w = SQLStorageBackend(st1_fname, mode='w')
    _store_via_container(be_1w, container, inp1, call_count=1,
                         backend_count=1)

    GeneralStorage._known_storages = {}
    storage = GeneralStorage(
        backend=SQLStorageBackend(st1_fname, mode='r'),
        class_info=_serialization,
        schema=_schema
    )
    container = _load_container(storage, sf)

    assert container.cv.func.call_count == 1
    os.remove(st1_fname)  # naughty user behavior

    with pytest.warns(UserWarning):
        assert container(inp1) == 'f'
    assert container.cv.func.call_count == 2
