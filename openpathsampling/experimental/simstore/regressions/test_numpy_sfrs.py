# This test checks a regression related to  storing results from storable
# functions that are of not basic types that SQLAlchemy knows. In
# particular, functions that returns NumPy arrays were storing the results
# without complaint, but could not re-load them as the correct type.
import pytest
import numpy as np

from openpathsampling.experimental.simstore import (
    StorableFunction, StorableFunctionResults, storable_function_find_uuids
)
from openpathsampling.experimental.simstore import SQLStorageBackend
from openpathsampling.experimental.simstore import GeneralStorage
from openpathsampling.experimental.simstore.class_info import (
    SerializationSchema, ClassInfo
)
from openpathsampling.netcdfplus import StorableObject
from openpathsampling.experimental.storage.ops_storage import (
    unsafe_ops_codecs, safe_ops_codecs
)
from openpathsampling.experimental.simstore.serialization_helpers import (
    default_find_uuids
)

from openpathsampling.experimental.simstore.test_utils import MockUUIDObject


@pytest.fixture
def minimal_sfr_serialization():
    schema = {
        'foo': [('normal_attr', 'int')],
        'storable_functions': [('json', 'json_obj')],
        'simulation_objects': [('json', 'json_obj')],
    }
    serialization = SerializationSchema(
        default_info=ClassInfo(
            table='simulation_objects',
            cls=StorableObject,
            serializer=unsafe_ops_codecs.simobj_serializer,
            deserializer=unsafe_ops_codecs.simobj_deserializer,
            safe_deserializer=safe_ops_codecs.simobj_serializer,
            find_uuids=default_find_uuids
        ),
        sfr_info=ClassInfo(
            table="function_results",
            cls=StorableFunctionResults,
            serializer=StorableFunctionResults.to_dict,
            deserializer=lambda x: x,  # deserializer not used
            find_uuids=default_find_uuids
        ),
        schema=schema,
        class_info_list=[
            ClassInfo(table='foo', cls=MockUUIDObject),
            ClassInfo(table='storable_functions',
                      cls=StorableFunction,
                      find_uuids=storable_function_find_uuids,
                      serializer=unsafe_ops_codecs.simobj_serializer,
                      deserializer=unsafe_ops_codecs.simobj_deserializer),
        ]
    )
    return schema, serialization

def stupid_func(inp):
    import numpy as np
    return np.array([1.0, 2.0, 3.0])

def test_numpy_sfrs(minimal_sfr_serialization, tmpdir):
    schema, serialization = minimal_sfr_serialization
    # create and save
    backend = SQLStorageBackend(tmpdir / 'tmp.db', mode='w')
    storage = GeneralStorage(backend, serialization, schema)

    inp = MockUUIDObject('foo', 1)
    func = StorableFunction(stupid_func).named('stupid')
    eval_result = func(inp)

    storage.save([func, inp])
    storage.close()

    # reload
    storage._known_storages = {}
    backend = SQLStorageBackend(tmpdir / 'tmp.db', mode='r')
    storage = GeneralStorage(backend, serialization, schema)

    sfs = list(storage.backend.table_iterator('storable_functions'))
    func = storage.load([sfs[0].uuid])[0]

    store_result = func(inp)

    assert store_result is not eval_result  # not same in memory
    assert isinstance(store_result, eval_result.__class__)
    np.testing.assert_array_equal(store_result, eval_result)


