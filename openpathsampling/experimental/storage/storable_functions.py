import collections
import warnings
import inspect
import types
from .tools import none_to_default
from .callable_codec import CallableCodec
from .serialization_helpers import get_uuid, has_uuid, default_find_uuids
from .class_info import ClassInfo
from .serialization_helpers import from_json_obj as deserialize_sim
from openpathsampling.netcdfplus import StorableNamedObject

import logging
logger = logging.getLogger(__name__)

class StorableFunctionResults(StorableNamedObject):
    """Cache of results from a function.

    The function will be associated with a cache like this, and this is the
    first place it will look for results (in order to avoid recalculating).
    """
    _sanity_checks = True
    def __init__(self, parent, parent_uuid):
        super(StorableFunctionResults, self).__init__()
        self.parent = parent
        self._parent_uuid = parent_uuid
        # TODO: someday, result_dict may need to be a cache that gets
        # emptied (thinking about memory concerns)
        self.result_dict = {}
        self.local_uuids = set([])

    @property
    def parent_uuid(self):
        if self.parent is not None:
            self._parent_uuid = get_uuid(self.parent)
        return self._parent_uuid

    def get_results_as_dict(self, uuid_items):
        """
        Parameters
        ----------
        items : Dict[str, obj]
            mapping of UUID to input item with that UUID

        Returns
        -------
        results : Dict[str, obj]
            mapping UUID to output for the input with that UUID
        missing : Set[str]
            set of UUIDs that could not be found (must be calculated)
        """
        uuids = set(uuid_items.keys())
        local_uuids = uuids & self.local_uuids
        missing_uuids = uuids - local_uuids

        results = {uuid: self.result_dict[uuid] for uuid in local_uuids}
        missing = {uuid: uuid_items[uuid] for uuid in missing_uuids}
        return results, missing

    def caching_sanity_check(self, uuid_results):
        for uuid, res in uuid_results.items():
            if uuid in self.local_uuids and self.result_dict[uuid] != res:
                # this only happens if we somehow get two different results
                # for the same input, which should not be possible
                # TODO: check which type of error this should be
                raise RuntimeError("Inconsistent results for storable"
                                   + "function.")

    def update(self, mapping):
        self.cache_results(mapping.result_dict)

    def cache_results(self, uuid_results):
        if self._sanity_checks:
            self.caching_sanity_check(uuid_results)
        self.result_dict.update(uuid_results)
        new_uuids = set(uuid_results.keys())
        self.local_uuids.update(new_uuids)

    def clear(self):
        self.result_dict.clear()
        self.local_uuids = set([])

    def __len__(self):
        return len(self.result_dict)

    def to_dict(self):
        return {
            'parent_uuid': self.parent_uuid,
            'result_dict': self.result_dict,
        }

    @classmethod
    def from_dict(cls, dct):
        obj = cls(dct['parent_uuid'])
        obj.result_dict = dct['result_dict']


class StorableFunction(StorableNamedObject):
    """Function wrapper, providing result caching and storage to disk.

    Parameters
    ----------
    func : Callable
    result_type
    store_store : Union[bool, None]
        Whether to store the source for this function. Default behavior
        (None) stores source for anything created in ``__main__`` that is
        not a lambda expression.
    """
    def __init__(self, func, result_type=None, store_source=None):
        super(StorableFunction, self).__init__()
        self.func = func
        self.source = None
        if store_source is None:
            is_lambda = (isinstance(func, types.FunctionType)
                         and func.__name__ == "<lambda>")
            store_source = func.__module__ == "__main__" and not is_lambda

        if store_source:
            try:
                self.source = inspect.getsource(func)
            except IOError:
                warnings.warn("Unable to get source for " + str(func))

        # self.input_type = input_type
        self.result_type = result_type
        self.local_cache = None  # set correctly by self.mode setter
        self._disk_cache = True
        self._handler = None
        self.mode = 'analysis'

    def set_handler(self, storage, override=False):
        if not override and self._handler is not None:
            raise RuntimeError("Handler for this StorableFunction has "
                               + "already been set:" + str(self._handler))
        self._handler = storage

    @property
    def has_handler(self):
        return self._handler is not None

    @property
    def disk_cache(self):
        return self._disk_cache

    @disk_cache.setter
    def disk_cache(self, value):
        self._disk_cache = value

    @property
    def mode(self):
        return self._mode

    @mode.setter
    def mode(self, value):
        allowed_values = ['no-caching', 'analysis', 'production']
        if value not in allowed_values:
            raise ValueError("Unknown mode: '%s'. Allowed options: %s" %
                             (value, allowed_values))

        if value == 'no-caching':
            self.local_cache = None
        else:
            self.local_cache = StorableFunctionResults(self, get_uuid(self))

        self._mode = value

    def to_dict(self):
        return {
            'func': self.func,  # made into JSON by CallableCodec
            'source': self.source,
            # 'input_type': self.input_type,
            'result_type': self.result_type,
        }

    @classmethod
    def from_dict(cls, dct):
        obj = cls(func=dct['func'],
                  # input_type=dct['input_type'],
                  result_type=dct['result_type'])
        if obj.source is None:
            obj.source = dct['source']  # may still be none

        return obj

    def preload_cache(self, storage=None):
        # TODO: add support for this to take a storage as argument
        if storage is None:
            storage = self._handler.storage

        uuid = get_uuid(self)
        cache_values = storage.backend.load_storable_function_table(uuid)
        self.local_cache.cache_results(cache_values)

    def is_singular(self, item):
        """Determine whether the input needs to be wrapped in a list.

        For performance (especially on analysis), all internal calculations
        are done assuming that the CV is called with a list of values. This
        method determines when the input was a list or a single item.

        Parameters
        ----------
        item : obj
            the input

        Returns
        -------
        bool :
            True if the item should be wrapped in a list
        """
        # the assumption here is that if you have a UUID, you're the type of
        # item we want. Need to override for snapshots-based CVs, where a
        # Trajectory has a UUID but is considered a list of Snapshots (maybe
        # we'll need to use isinstance(item, Snapshot) instead?)
        return has_uuid(item)

    def _eval(self, uuid_items):
        if self.func is None and uuid_items:
            raise RuntimeError("No function attached to %s. Can not "
                               + "evaluate for %s." % (self, uuid_items))

        results = {uuid: self.func(item)
                   for uuid, item in uuid_items.items()}
        return results, {}

    def _get_cached(self, uuid_items):
        return self.local_cache.get_results_as_dict(uuid_items)

    def _get_storage(self, uuid_items):
        if not self._handler:
            return {}, uuid_items

        return self._handler.get_function_results(get_uuid(self), uuid_items)

    def __call__(self, items):
        # important: implementation is that we always try to take an
        # iterable of items ... this makes, e.g., applying a CV to a
        # trajectory (especially in analysis) MUCH faster
        singular = False
        if self.is_singular(items):
            items = [items]
            singular = True

        uuid_items = {get_uuid(item): item for item in items}
        # TODO: add preprocessing here? if needed?

        cache_mode_order = {  # tuples of (function, add_to_cache)
            'analysis': [(self._get_cached, False),
                         (self._get_storage, True),
                         (self._eval, True)],
            'production': [(self._get_cached, False),
                           (self._eval, True)],
            'no-caching': [(self._eval, False)]
        }[self.mode]

        missing = uuid_items
        result_dict = {}
        for stage, do_caching in cache_mode_order:
            logger.debug(stage, missing)
            stage_results, missing = stage(missing)
            result_dict.update(stage_results)
            if do_caching:
                self.local_cache.cache_results(stage_results)

            if not missing:
                break

        return_list = [result_dict[get_uuid(item)] for item in items]
        # result = self.postprocess(return_list)  # TODO if needed?
        result = return_list

        if singular:
            result = result[0]

        return result


def storable_function_find_uuids(obj, cache_list):
    uuids, new_objects = default_find_uuids(obj, cache_list)
    func_results = obj.local_cache
    uuids.update({get_uuid(func_results): func_results})
    return uuids, new_objects


class SFRClassInfo(ClassInfo):
    def __init__(self, table='function_results',
                 cls=StorableFunctionResults,
                 serializer=StorableFunctionResults.to_dict,
                 deserializer=None, find_uuids=default_find_uuids):
        deserializer = none_to_default(deserializer, lambda x: x)
        super(self, SFRClassInfo).__init__(table=table, cls=cls,
                                           serializer=serializer,
                                           deserializer=deserlializer,
                                           find_uuids=find_uuids)


class StorageFunctionHandler(object):
    def __init__(self, storage, functions=None, other_codecs=None,
                 codec_settings=None):
        self.storage = storage
        self.codecs = none_to_default(other_codecs, [])
        functions = none_to_default(functions, [])
        codec_settings = none_to_default(codec_settings, {})
        self.callable_codec = None
        self._codec_settings = None
        self.codec_settings = codec_settings  # sets callable_codec
        self.canonical_functions = {}
        self.all_functions = collections.defaultdict(list)
        for func in functions:
            self.register_function(func)

    @property
    def codec_settings(self):
        return self._codec_settings

    @codec_settings.setter
    def codec_settings(self, settings):
        if self._codec_settings != settings:
            self._codec_settings = settings
            self.callable_codec = CallableCodec(settings)

    def _make_classinfo(self, func):
        pass

    def register_function(self, func, add_table=True):
        func_uuid = get_uuid(func)
        # register as a canonical function
        if func_uuid not in self.canonical_functions:
            logger.debug("Registering new function: %s" % func_uuid)
            self.canonical_functions[func_uuid] = func
            if add_table:
                self.storage.backend.register_storable_function(
                    table_name=func_uuid,
                    result_type=func.result_type
                )

        # register will all_functions
        is_registered = any([func is registered
                             for family in self.all_functions.values()
                             for registered in family])
        if not is_registered:
            self.all_functions[func_uuid].append(func)

        # set handler
        if not func.has_handler:
            func.set_handler(self)

    def clear_non_canonical(self):
        self.all_functions = collections.defaultdict(list)
        for func in self.functions:
            # re-registering the canonical func will put it in all_functions
            self.register_function(func)

    @property
    def functions(self):
        return list(self.canonical_functions.values())

    def update_cache(self, function_results):
        self.register_function(function_results.parent)
        func = self.canonical_functions[function_results.parent_uuid]
        func.local_cache.update(function_results)

    def get_function_results(self, func_uuid, uuid_items):
        backend = self.storage.backend
        uuids = set(uuid_items.keys())
        uuid_map = backend.load_storable_function_results(func_uuid, uuids)
        missing = uuids - set(uuid_map.keys())
        missing_map = {uuid: uuid_items[uuid] for uuid in missing}
        return uuid_map, missing_map

