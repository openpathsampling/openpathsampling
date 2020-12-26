import functools

import dask
import distributed
from distributed.protocol.serialize import register_serialization_family

import cloudpickle as cp

from openpathsampling.experimental.simstore.memory_backend import \
        MemoryStorageBackend
from openpathsampling.experimental.simstore.sql_backend import \
        SQLStorageBackend
from openpathsampling.experimental.storage.ops_storage import Storage
from openpathsampling.experimental.simstore.serialization_helpers import (
    has_uuid, get_uuid
)

from openpathsampling.experimental.simstore.tools import none_to_default

from openpathsampling.netcdfplus import StorableNamedObject

import os

def ops_dumps(obj):
    header = {'serializer': 'openpathsampling'}
    if not has_uuid(obj):
        # quick exit
        raise NotImplementedError()

    storage = Storage.from_backend(MemoryStorageBackend())
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
    storage = Storage.from_backend(backend)
    for func in storage.storable_functions:
        func.preload_cache(storage)
    obj = storage.load([dct['uuid']])[0]
    return obj


class OPSTask(StorableNamedObject):
    """Simple wrapper for task, args, and kwargs to make storable"""
    def __init__(self, task, args, kwargs):
        super().__init__()
        self.task = serialize_task(task)
        self.args = args
        self.kwargs = kwargs

    @classmethod
    def from_dict(cls, dct):
        obj = cls.__new__(cls)
        obj.task = dct['task']
        obj.args = dct['args']
        obj.kwargs = dct['kwargs']
        return obj

    def call_with_args(self, *args, **kwargs):
        task = deserialize_task(self.task)
        return task(*args, **kwargs)


    def __call__(self):
        task = deserialize_task(self.task)
        return task(*self.args, **self.kwargs)


class PicklableOPSTask(object):
    def __init__(self, task, args=None, kwargs=None):
        args = none_to_default(args, [])
        kwargs = none_to_default(kwargs, {})
        ops_task = OPSTask(task, args, kwargs)
        self.backend = MemoryStorageBackend()
        storage = Storage.from_backend(self.backend)
        self.task_uuid = get_uuid(ops_task)
        storage.save(ops_task)

    @classmethod
    def reconstruct(cls, task_uuid, backend):
        obj = cls.__new__(cls)
        obj.task_uuid = task_uuid
        obj.backend = backend
        obj.backend.mode = 'r'
        return obj

    def __reduce__(self):
        return (PicklableOPSTask.reconstruct, (self.task_uuid, self.backend))

    def _reconstruct_task(self):
        storage = Storage.from_backend(self.backend)
        ops_task = storage.load([self.task_uuid])[0]
        return ops_task

    def call_with_args(self, *args, **kwargs):
        ops_task = self._reconstruct_task()
        # print(f"call_with_args: {ops_task}: args={args} kwargs={kwargs}")
        return ops_task.call_with_args(*args, **kwargs)

    def __call__(self):
        ops_task = self._reconstruct_task()
        return ops_task()


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


class SerialScheduler(object):
    def store_results(self, storage, results):
        storage.save(results)

    def wrap_task_only(self, task):
        return task

    def wrap_task(self, task):
        return task

    def finalize(self):
        pass


class DaskDistributedScheduler(object):
    def __init__(self, client=None):
        if client is None:
            client = distributed.Client(
                serializers=['openpathsampling', 'dask', 'pickle',  ],
                deserializers=['openpathsampling', 'dask', 'pickle', ],
            )
        self.client = client
        self._final_future = None

    def store_results(self, storage, result_future):
        """Store a result

        Parameters
        ----------
        storage : :class:`.Storage`
            OPS SimStore storage to save results in
        result_future : dask.Future
            Result to store. Note that this MUST be a dask Future! Regular
            OPS objects will not retain UUIDs if pickled.
        """
        # TODO: make this work with both futures and regular OPS objects?
        # print(f"About to save {result_future} to {storage}")
        ops_task = PicklableOPSTask(storage.save)
        fut = self.client.submit(ops_task.call_with_args, result_future)
        distributed.fire_and_forget(fut)
        self._final_future = fut

    def wrap_task_only(self, task):
        @functools.wraps(task)
        def inner(*args, **kwargs):
            ops_task = PicklableOPSTask(task)
            d_args= [dask.delayed(arg) for arg in args]
            d_kwargs = {key: dask.delayed(kwarg)
                        for key, kwarg in kwargs.items()}
            result = self.client.submit(ops_task.call_with_args, *args,
                                        **kwargs, pure=False)
            self._final_future = result
            return result
        return inner

    def wrap_task(self, task):
        @functools.wraps(task)
        def inner(*args, **kwargs):
            ops_task = PicklableOPSTask(task, args, kwargs)
            # technically, don't think we need pure=False, since every
            # ops_task will be a different object in memory
            result = self.client.submit(ops_task, pure=False)
            self._final_future = result
            return result

        return inner

    def finalize(self):
        # ensure that we run the last future
        if self._final_future is not None:
            _ = self._final_future.result()
