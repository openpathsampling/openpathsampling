import collections
import warnings
import inspect
import types

import numpy as np

from .tools import none_to_default
from .callable_codec import CallableCodec
from .serialization_helpers import get_uuid, has_uuid, default_find_uuids
from .class_info import ClassInfo
from .serialization_helpers import from_json_obj as deserialize_sim
from openpathsampling.netcdfplus import StorableNamedObject

import logging
logger = logging.getLogger(__name__)


class Processor(StorableNamedObject):
    """Storable function pre/post processors"""
    def __init__(self, name, func, stage):
        super(Processor, self).__init__()
        self.name = name
        self.func = func
        stages = ['item-pre', 'list-pre', 'item-post', 'list-post']
        if stage in stages:
            self.stage = stage
        else:
            raise ValueError("Unknown state: '%s'. Allowed stages are: %s"
                             % (stage, stages))

    def __call__(self, inputs):
        return self.func(inputs)

requires_lists_pre = Processor(name='requires_lists_pre',
                               func=lambda values: [list(values)],
                               stage='list-pre')

requires_lists_post = Processor(name='requires_lists_post',
                                func=lambda results: results[0],
                                stage='list-post')

def _scalarize_singletons(values):
    """Post-processing to scalarize singletons within a list of results.

    If each snapshot returns a list, this converts length-1 lists into
    scalars.

    .. note::
    The current implementation only works on NumPy arrays, and will ignore
    all other data.
    """
    if isinstance(values, np.ndarray):
        shape = tuple(n for n in values.shape if n != 1)
        if shape == tuple():
            values = values.__float__()
        else:
            values.shape = shape

        # shape = values.shape
        # if len(shape) > 1 and shape[1] == 1:
            # new_shape = tuple([shape[0]] + list(shape)[2:])
            # values.shape = new_shape
    return values

scalarize_singletons = Processor(name='scalarize_singletons',
                                 func=_scalarize_singletons,
                                 stage='item-post')

wrap_numpy = Processor(name="wrap_numpy",
                       func=lambda values: np.array(values),
                       stage='list-post')


class StorableFunctionConfig(StorableNamedObject):
    """Manages pre/post processing options for CVs.

    This allows simple :class:`.Processor` instances to be combined in order
    to facilitate more complicated pre- and post-processing of inputs to the
    evaluation stage.
    """
    def __init__(self, processors=None):
        if processors is None:
            processors = []
        self.processors = []
        self.processor_dict = {}
        self.list_preprocessors = []
        self.item_preprocessors = []
        self.list_postprocessors = []
        self.item_postprocessors = []
        for proc in processors:
            self.register(proc)

    def _build_processor_lists(self):
        stage_to_list = {'item-pre': self.item_preprocessors,
                         'item-post': self.item_postprocessors,
                         'list-pre': self.list_preprocessors,
                         'list-post': self.list_postprocessors}
        for proc_list in stage_to_list.values():
            del proc_list[:]  # https://stackoverflow.com/a/1400622

        for proc in self.processors:
            stage_to_list[proc.stage].append(proc)

    def register(self, processor):
        if processor.name in self.processor_dict:
            self.deregister(processor.name)
        self.processor_dict[processor.name] = processor
        self.processors.append(processor)
        self._build_processor_lists()

    def deregister(self, processor, error_if_missing=True):
        name = processor if isinstance(processor, str) else processor.name
        to_remove = self.processor_dict.pop(name, None)
        if to_remove is None:
            if error_if_missing:
                raise KeyError(processor)
            else:
                return
        self.processors.remove(to_remove)
        self._build_processor_lists()

    @staticmethod
    def _process(processors, inputs):
        for proc in processors:
            inputs = proc(inputs)
        return inputs

    def item_preprocess(self, values):
        return self._process(self.item_preprocessors, values)

    def list_preprocess(self, values):
        return self._process(self.list_preprocessors, values)

    def item_postprocess(self, values):
        return self._process(self.item_postprocessors, values)

    def list_postprocess(self, values):
        return self._process(self.list_postprocessors, values)


class StorableFunctionResults(StorableNamedObject):
    """Cache of results from a function.

    The function will be associated with a cache like this, and this is the
    first place it will look for results (in order to avoid recalculating).
    """
    _sanity_checks = False  # off by default; can't work with ndarray
    def __init__(self, parent, parent_uuid):
        # TODO: why does this require parent_uuid still?
        super(StorableFunctionResults, self).__init__()
        self.parent = parent
        self._parent_uuid = parent_uuid
        # TODO: someday, result_dict may need to be a cache that gets
        # emptied or an LRU cache (thinking about memory concerns)
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


def _storage_exception_msg(func, storage, exc):
    func_id = str(get_uuid(func))
    if func.is_named:
        func_id += " (" + func.name + ")"
    try:
        storage_id = " in " + str(storage.identifier) + "."
    except:
        storage_id = "."

    msg = ("A problem occured while loading disk-cached results. "
           "Function " + func_id + storage_id + "\n" + str(type(exc))
           + str(exc))
    return msg


class StorableFunction(StorableNamedObject):
    """Function wrapper, providing result caching and storage to disk.

    Parameters
    ----------
    func : Callable
    store_store : Union[bool, None]
        Whether to store the source for this function. Default behavior
        (None) stores source for anything created in ``__main__`` that is
        not a lambda expression.
    """
    def __init__(self, func, func_config=None, store_source=None,
                 period_min=None, period_max=None,**kwargs):
        super(StorableFunction, self).__init__()
        self.func = func
        self.source = None
        self.kwargs = kwargs
        if func_config is None:
            func_config = StorableFunctionConfig()
        self.func_config = func_config
        if store_source is None:
            is_lambda = (isinstance(func, types.FunctionType)
                         and func.__name__ == "<lambda>")
            store_source = func.__module__ == "__main__" and not is_lambda

        if store_source:
            try:
                self.source = inspect.getsource(func)
            except IOError:
                warnings.warn("Unable to get source for " + str(func))

        self.local_cache = None  # set correctly by self.mode setter
        self._disk_cache = True
        self._handlers = set([])
        self._check_period(period_min, period_max)
        self.period_min = period_min
        self.period_max = period_max
        self._modes = {  # tuples of (function, add_to_cache)
            'analysis': [(self._get_cached, False),
                         (self._get_storage, True),
                         (self._eval, True)],
            'production': [(self._get_cached, False),
                           (self._eval, True)],
            'no-caching': [(self._eval, False)]
        }
        self.mode = 'analysis'

    @staticmethod
    def _check_period(period_min, period_max):
        is_not_periodic = period_min is None and period_max is None
        is_periodic = period_min is not None and period_max is not None
        is_error = not (is_periodic or is_not_periodic)
        if is_error:
            raise ValueError("Periodic functions must have upper and "
                             "lower bounds. This function has period_min "
                             + str(period_min) + " and period_max "
                             + str(period_max) + ".")
        return is_periodic

    @property
    def is_periodic(self):
        return self._check_period(self.period_min, self.period_max)


    def add_handler(self, storage, override=False):
        self._handlers.add(storage)

    def remove_handler(self, handler):
        self._handlers.discard(handler)

    @property
    def has_handler(self):
        return len(self._handlers) > 0

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
        allowed_values = list(self._modes)
        if value not in allowed_values:
            raise ValueError("Unknown mode: '%s'. Allowed options: %s" %
                             (value, allowed_values))

        if value == 'no-caching':
            self.local_cache = None
            self.disk_cache = False
        else:
            self.local_cache = StorableFunctionResults(self, get_uuid(self))

        self._mode = value

    def to_dict(self):
        return {
            'func': self.func,  # made into JSON by CallableCodec
            'source': self.source,
            'kwargs': self.kwargs,
        }

    @classmethod
    def from_dict(cls, dct):
        source = dct.pop('source')
        kwargs = dct.pop('kwargs')
        dct['store_source'] = False
        obj = cls(**dct, **kwargs)
        if obj.source is None:
            obj.source = source  # may still be none

        return obj

    def preload_cache(self, storage=None):
        if storage is None:
            storages = [h.storage for h in self._handlers]
        else:
            storages = [storage]

        uuid = get_uuid(self)
        for storage in storages:
            cache = storage.backend.load_storable_function_table(uuid)
            self.local_cache.cache_results(cache)

    def is_scalar(self, item):
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
        values = list(uuid_items.values())
        preprocessed = self.func_config.item_preprocess(values)
        preprocessed = self.func_config.list_preprocess(preprocessed)
        values = [self.func(item, **self.kwargs) for item in preprocessed]
        postprocessed = [self.func_config.item_postprocess(val)
                         for val in values]
        results = dict(zip(uuid_items.keys(), postprocessed))
        return results, {}

    def _get_cached(self, uuid_items):
        return self.local_cache.get_results_as_dict(uuid_items)

    def _get_storage(self, uuid_items):
        if not self.has_handler:
            return {}, uuid_items

        missing = uuid_items
        found = {}
        my_uuid = get_uuid(self)
        for handler in self._handlers:
            try:
                uuid_map, missing = handler.get_function_results(
                    my_uuid, missing
                )
            except Exception as e:
                # for any error, we just warn -- can be correctly calculated
                # in eval
                msg = _storage_exception_msg(self, handler.storage, e)
                warnings.warn(msg)
                uuid_map = {}

            found.update(uuid_map)

        return found, missing

    def __call__(self, items):
        # important: implementation is that we always try to take an
        # iterable of items ... this makes, e.g., applying a CV to a
        # trajectory (especially in analysis) MUCH faster
        scalar = False
        if self.is_scalar(items):
            items = [items]
            scalar = True

        uuid_items = {get_uuid(item): item for item in items}
        # TODO: add preprocessing here? if needed?

        cache_mode_order = self._modes[self.mode]

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

        result = [result_dict[get_uuid(item)] for item in items]

        result = self.func_config.list_postprocess(result)

        if scalar:
            # TODO: watch to see if this try-catch causes performance issues
            try:
                result = result[0]
            except IndexError:
                # this means that the result is already a scalar
                pass

        return result


def storable_function_find_uuids(obj, cache_list):
    # TODO: it should be possible to remove this at some point
    uuids, new_objects = default_find_uuids(obj, cache_list)
    if obj.local_cache is not None:
        uuids[get_uuid(obj.local_cache)] = obj.local_cache

    return uuids, new_objects


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

    def register_function(self, func, example_result=None):
        func_uuid = get_uuid(func)

        # add table to backend if needed
        needs_table = not self.storage.backend.has_table(func_uuid)
        add_table = needs_table and (example_result is not None)
        if needs_table and not add_table:
            logger.info("Result type unknown; unable to create table")
        elif add_table:
            identify = self.storage.type_identification.identify
            result_type = identify(example_result)
            self.storage.backend.register_storable_function(
                table_name=func_uuid,
                result_type=result_type
            )

        # register as a canonical function
        if func_uuid not in self.canonical_functions:
            logger.debug("Registering new function: %s" % func_uuid)
            self.canonical_functions[func_uuid] = func

        # register with all_functions
        is_registered = any([func is registered
                             for family in self.all_functions.values()
                             for registered in family])
        if not is_registered:
            self.all_functions[func_uuid].append(func)

        if add_table or not needs_table:
            func.add_handler(self)

    def close(self):
        for uuid, copies in self.all_functions.items():
            for func in copies:
                func.remove_handler(self)

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
