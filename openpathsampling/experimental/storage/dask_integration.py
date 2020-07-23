import distributed
from distributed.protocol.serialize import register_serialization_family

import cloudpickle as cp

from openpathsampling.experimental.storage.memory_backend import \
        MemoryStorageBackend
from openpathsampling.experimental.storage.sql_backend import \
        SQLStorageBackend
from openpathsampling.experimental.storage.ops_storage import OPSStorage
from openpathsampling.experimental.storage.serialization_helpers import (
    has_uuid, get_uuid
)

import os

def ops_dumps(obj):
    header = {'serializer': 'openpathsampling'}
    if not has_uuid(obj):
        # quick exit
        raise NotImplementedError()

    storage = OPSStorage.from_backend(MemoryStorageBackend())
    storage.save(obj)
    uuid = get_uuid(obj)
    frames = [cp.dumps({'uuid': uuid, 'backend': storage.backend})]
    return header, frames

def ops_loads(header, frames):
    if len(frames) > 1:
        frame = ''.join(frames)
    else:
        frame = frames[0]
    dct = cp.loads(frame)
    backend = dct['backend']
    backend.mode = 'r'
    storage = OPSStorage.from_backend(backend)
    for func in storage.storable_functions:
        func.preload_cache(storage)
    obj = storage.load([dct['uuid']])[0]
    return obj


def serialize_task(task):
    # allows serialization of: (1) UUID objects; (2) methods of UUID object;
    # (3) bare functions
    if has_uuid(task):
        func = None
        obj = task
    elif hasattr(task, '__self__') and has_uuid(task.__self__):
        func = task.__name__
        obj = task.__self__
    else:
        func = task
        obj = None
    return obj, func

def deserialize_task(ser_task):
    obj, func = ser_task
    if obj is None:
        task = func
    elif func is None:
        task = obj
    else:
        task = getattr(obj, func)
    return task

register_serialization_family('openpathsampling', ops_dumps, ops_loads)


def _save_results(storage, results):
    print(f"Using storage {storage} to save {results}")
    storage.save(results)

def _remote_task(serialized_task, *args, **kwargs):
    task = deserialize_task(serialized_task)
    result = task(*args, **kwargs)
    return result


class SerialScheduler(object):
    def store_results(self, filename, results):
        _save_results(storage, results)

    def wrap_task(task):
        return task


class DaskDistributedScheduler(object):
    def __init__(self, client=None):
        if client is None:
            client = distributed.Client(
                serializers=['openpathsampling', 'dask', 'pickle',  ],
                deserializers=['openpathsampling', 'dask', 'pickle', ],
            )
        self.client = client
        self._final_future = None

    def store_results(self, storage, results):
        print(f"About to save {results} to {storage}")
        fut = self.client.submit(_save_results, storage, results)
        distributed.fire_and_forget(fut)
        self._final_future = fut

    def wrap_task(self, task):
        def inner(*args, **kwargs):
            ser_task = serialize_task(task)
            result = self.client.submit(_remote_task,
                                        ser_task,
                                        pure=False,
                                        *args,
                                        **kwargs)
            return result

        return inner

    def finalize(self):
        # ensure that we run the last future
        if self._final_future is not None:
            _ = self._final_future.result()
