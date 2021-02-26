"""
A simple storage interface for simulation objects and data objects.

Reserved words
--------------

Table names:

* ``uuid``
* ``tables``
* ``metadata``

Column names:

* ``uuid``
* ``idx``
"""

import logging
import collections
from collections import abc
import itertools

from . import tools
from .serialization_helpers import get_uuid, get_all_uuids
from .serialization_helpers import get_all_uuids_loading
from .serialization_helpers import get_reload_order
# from .serialization import Serialization
from .proxy import ProxyObjectFactory, GenericLazyLoader
from .storable_functions import StorageFunctionHandler, StorableFunction
from .tags_table import TagsTable
from .type_ident import STANDARD_TYPING

try:
    basestring
except NameError:
    basestring = str

logger = logging.getLogger(__name__)

# these two tables are required in *all* schema
universal_schema = {
    'uuid': [('uuid', 'uuid'), ('table', 'int'), ('row', 'int')],
    'tables': [('name', 'str'), ('idx', 'int'), ('module', 'str'),
               ('class_name', 'str')],
    'tags': [('name', 'str'), ('content', 'uuid')]
}


from openpathsampling.netcdfplus import StorableNamedObject
class GeneralStorage(StorableNamedObject):
    _known_storages = {}
    def __init__(self, backend, class_info, schema=None,
                 simulation_classes=None, fallbacks=None, safemode=False):
        super().__init__()
        GeneralStorage._known_storages[backend.identifier] = self
        self.backend = backend
        self.schema = schema.copy()
        self.class_info = class_info.copy()
        self._safemode = None
        self.safemode = safemode
        self._sf_handler = StorageFunctionHandler(storage=self)
        self.type_identification = STANDARD_TYPING  # TODO: copy
        # TODO: implement fallback
        self.fallbacks = tools.none_to_default(fallbacks, [])

        self.simulation_classes = tools.none_to_default(simulation_classes,
                                                        {})

        # self._pseudo_tables = {table_name: dict()
                               # for table_name in self.simulation_classes}
        self._simulation_objects = {}
        self._pseudo_tables = {table_name: PseudoTable()
                               for table_name in self.simulation_classes}
        self._pseudo_tables['misc_simulation'] = PseudoTable()

        self._storage_tables = {}  # stores .steps, .snapshots
        # self.serialization = Serialization(self)
        self.proxy_factory = ProxyObjectFactory(self, self.class_info)
        self.cache = MixedCache()
        if self.schema is None:
            self.schema = backend.schema
        self.cache = MixedCache({})  # initial empty cache so it exists
        self.initialize_with_mode(self.mode)
        self.tags = TagsTable(self)
        self._simulation_objects = self._cache_simulation_objects()
        self.cache = MixedCache(self._simulation_objects)
        self._stashed = []
        self._reset_fixed_cache()

    @property
    def mode(self):
        return self.backend.mode

    @property
    def safemode(self):
        return self._safemode

    @safemode.setter
    def safemode(self, value):
        if value is self._safemode:
            return
        self.class_info.set_safemode(value)

    def initialize_with_mode(self, mode):
        if mode == 'r' or mode == 'a':
            self.register_schema(self.schema, class_info_list=[],
                                 read_mode=True)
            missing = {k: v for k, v in self.backend.schema.items()
                       if k not in self.schema and k not in universal_schema}
            self.schema.update(missing)
            table_to_class = self.backend.table_to_class
            self._load_missing_info_tables(table_to_class)
            sim_objs = {get_uuid(obj): obj
                        for obj in self.simulation_objects}
            sim_objs.update({get_uuid(obj): obj
                             for obj in self.storable_functions})
            self._simulation_objects.update(sim_objs)
            self._update_pseudo_tables(sim_objs)

        elif mode == 'w':
            self.register_schema(self.schema, class_info_list=[])

    def _load_missing_info_tables(self, table_to_class):
        missing_info_tables = [tbl for tbl in self.schema
                               if tbl not in self.class_info.tables]
        n_missing = len(missing_info_tables)
        logger.info("Missing info from %d dynamically-registered tables",
                    n_missing)
        classes = [table_to_class[tbl] for tbl in missing_info_tables]
        self.register_from_tables(missing_info_tables, classes)
        missing_info_tables = [tbl for tbl in self.schema
                               if tbl not in self.class_info.tables]
        logger.info("Successfully registered %d missing tables",
                    n_missing - len(missing_info_tables))

        if missing_info_tables:
            raise RuntimeError("Unable to register existing database "
                               + "tables: " + str(missing_info_tables))

    def register_from_tables(self, table_names, classes):
        # override in subclass to handle special lookups
        pass

    def stash(self, objects):
        objects = tools.listify(objects)
        self._stashed.extend(objects)

    def close(self):
        # TODO: should sync on close
        self.backend.close()
        self._sf_handler.close()
        for fallback in self.fallbacks:
            fallback.close()

    def register_schema(self, schema, class_info_list,
                        backend_metadata=None, read_mode=False):
        # check validity
        self.class_info.register_info(class_info_list, schema)
        # for info in class_info_list:
            # info.set_defaults(schema)
            # self.class_info.add_class_info(info)

        schema_types = [type_str for attr_list in schema.values()
                        for _, type_str in attr_list]
        for type_str in schema_types:
            backend_type = self.class_info.backend_type(type_str)
            self.backend.register_type(type_str, backend_type)

        if not read_mode or self.backend.table_to_class == {}:
            table_to_class = {table: self.class_info[table].cls
                              for table in schema
                              if table not in ['uuid', 'tables']}
            # here's where we add the class_info to the backend
            self.backend.register_schema(schema, table_to_class,
                                         backend_metadata)

        self.schema.update(schema)
        for table in self.schema:
            self._storage_tables[table] = StorageTable(self, table)
        # self.serialization.register_serialization(schema, self.class_info)

    def register_from_instance(self, lookup, obj):
        raise NotImplementedError("No way to register from an instance")

    def register_missing_tables_for_objects(self, uuid_obj_dict):
        # missing items are handled by the special_lookup
        lookup_examples = set([])
        for obj in uuid_obj_dict.values():
            lookup = self.class_info.lookup_key(obj)
            if lookup not in lookup_examples:
                self.register_from_instance(lookup, obj)
                lookup_examples |= {lookup}

    def filter_existing_uuids(self, uuid_dict):
        existing = self.backend.load_uuids_table(
            uuids=list(uuid_dict.keys()),
            ignore_missing=True
        )

        # "special" here indicates that we always try to re-save these, even
        # if they've already been saved once. This is (currently) necessary
        # for objects that contain mutable components (such as the way
        # storable functions contain their results)
        # TODO: make `special` customizable
        special = set(self._sf_handler.canonical_functions.keys())

        for uuid_row in existing:
            uuid = uuid_row.uuid
            if uuid not in special:
                del uuid_dict[uuid_row.uuid]

        return uuid_dict


    def _uuids_by_table(self, input_uuids, cache, get_table_name):
        # find all UUIDs we need to save with this object
        logger.debug("Listing all objects to save")
        uuids = {}
        for uuid, obj in input_uuids.items():
            uuids.update(get_all_uuids(obj, known_uuids=cache,
                                       class_info=self.class_info))

        logger.debug("Found %d objects" % len(uuids))
        logger.debug("Deproxying proxy objects")
        uuids = self._unproxy_lazies(uuids)

        logger.debug("Checking if objects already exist in database")
        uuids = self.filter_existing_uuids(uuids)

        # group by table, then save appropriately
        # by_table; convert a dict of {uuid: obj} to {table: {uuid: obj}}
        by_table = collections.defaultdict(dict)
        by_table.update(tools.dict_group_by(uuids,
                                            key_extract=get_table_name))
        return by_table

    def _unproxy_lazies(self, uuid_mapping):
        lazies = [uuid for uuid, obj in uuid_mapping.items()
                  if isinstance(obj, GenericLazyLoader)]
        logger.debug("Found " + str(len(lazies)) + " objects to deproxy")
        loaded = self.load(lazies, allow_lazy=False)
        uuid_mapping.update({get_uuid(obj): obj for obj in loaded})
        return uuid_mapping


    def save(self, obj_list, use_cache=True):
        if type(obj_list) is not list:
            obj_list = [obj_list]

        cache = self.cache if use_cache else {}
        # TODO: convert the whole .save process to something based on the
        # class_info.serialize method (enabling per-class approaches for
        # finding UUIDs, which will be a massive serialization speed-up
        # self.class_info.serialize(obj, storage=self)
        # check if obj is in DB (maybe this can be removed?)
        logger.debug("Starting save")
        input_uuids = {get_uuid(obj): obj for obj in obj_list}
        input_uuids = self.filter_existing_uuids(input_uuids)
        if not input_uuids:
            return  # exit early if everything is already in storage

        # check default table for things to register; register them
        # TODO: move to function: self.register_missing(by_table)
        # TODO: convert to while?
        get_table_name = lambda uuid, obj_: self.class_info[obj_].table
        by_table = self._uuids_by_table(input_uuids, cache, get_table_name)
        old_missing = {}
        while '__missing__' in by_table:
            # __missing__ is a special result returned by the
            # ClassInfoContainer if this is object is expected to have a
            # table, but the table doesn't exist (e.g., for dynamically
            # added tables)
            missing = by_table.pop('__missing__')
            if missing == old_missing:
                raise RuntimeError("Unable to register: " + str(missing))
            missing = self._unproxy_lazies(missing)
            logger.info("Registering tables for %d missing objects",
                        len(missing))
            self.register_missing_tables_for_objects(missing)
            missing_by_table = tools.dict_group_by(missing, get_table_name)
            logger.info("Registered %d new tables: %s",
                        len(missing_by_table),
                        str(list(missing_by_table.keys())))
            by_table.update(missing_by_table)
            # search for objects inside the objects we just registered
            next_by_table = self._uuids_by_table(missing, cache,
                                                 get_table_name)
            for table, uuid_dict in next_by_table.items():
                by_table[table].update(uuid_dict)
            old_missing = missing

        # TODO: move to function self.store_sfr_results(by_table)
        self.save_function_results()  # always for canonical
        has_sfr = (self.class_info.sfr_info is not None
                   and self.class_info.sfr_info.table in by_table)
        if has_sfr:
            func_results = by_table.pop(self.class_info.sfr_info.table)
            logger.info("Saving results from %d storable functions",
                        len(func_results))
            for result in func_results.values():
                func = result.parent
                self.save_function_results(func)

        # TODO: add simulation objects to the cache

        # this is the actual serialization
        logger.debug("Filling %d tables: %s", len(by_table),
                     str(list(by_table.keys())))
        for table in by_table:
            logger.debug("Storing %d objects to table %s",
                         len(by_table[table]), table)
            serialize = self.class_info[table].serializer
            storables_list = [serialize(o) for o in by_table[table].values()]
            self.backend.add_to_table(table, storables_list)
            # special handling for simulation objects
            if table == 'simulation_objects':
                self._update_pseudo_tables(by_table[table])
                self._simulation_objects.update(by_table[table])
                self._reset_fixed_cache()
            logger.debug("Storing complete")

    def save_function_results(self, funcs=None):
        # TODO: move this to sf_handler; where the equivalent load happens

        # no equivalent load because user has no need -- any loading can be
        # done by func, either as func(obj) or func.preload_cache()
        if funcs is None:
            funcs = list(self._sf_handler.canonical_functions.values())

        funcs = tools.listify(funcs)
        for func in funcs:
            # TODO: XXX This is where we need to use type identification
            # 1. check if the associated function is registered already
            # 2. if not, extract type from the first function value
            self._sf_handler.update_cache(func.local_cache)
            result_dict = func.local_cache.result_dict
            table_name = get_uuid(func)
            if not self.backend.has_table(table_name):
                if result_dict:
                    example = next(iter(result_dict.values()))
                else:
                    example = None
                self._sf_handler.register_function(func,
                                                   example_result=example)

            if result_dict:
                self.backend.add_storable_function_results(
                    table_name=table_name,
                    result_dict=result_dict
                )
            self._reset_fixed_cache()

    def load(self, input_uuids, allow_lazy=True, force=False):
        """
        Parameters
        ----------
        input_uuids : List[str]
        allow_lazy : bool
            whether to allow lazy proxy objects
        force : bool
            force reloading this object even if it is already cached (used
            for deproxying a lazy proxy object)
        """
        # loading happens in 4 parts:
        # 1. Get UUIDs that need to be loaded
        # 2. Make lazy-loading proxy objects
        # 3. Identify the order in which we deserialize
        # 4. Deserialize
        # set force=True to make it reload this full object (used for
        # loading a lazy-loaded object)
        if isinstance(input_uuids, basestring):
            # TEMP: remove; for now, prevents my stupidity
            raise RuntimeError("David, you forgot to wrap UUID in list")

        logger.debug("Starting to load %d objects", len(input_uuids))
        if force:
            self.cache.delete_items(input_uuids)

        results = {uuid: self.cache[uuid] for uuid in input_uuids
                   if uuid in self.cache}
        uuid_list = [uuid for uuid in input_uuids if uuid not in self.cache]
        logger.debug("Getting internal structure of %d non-cached objects",
                     len(uuid_list))
        to_load, lazy_uuids, dependencies, uuid_to_table = \
                get_all_uuids_loading(uuid_list=uuid_list,
                                      backend=self.backend,
                                      schema=self.schema,
                                      existing_uuids=self.cache,
                                      allow_lazy=allow_lazy)
        logger.debug("Loading %d objects; creating %d lazy proxies",
                     len(to_load), len(lazy_uuids))

        # to_load : List (table rows from backend)
        # lazy : Set[str] (lazy obj UUIDs)
        # dependencies : Dict[str, List[str]] (map UUID to contained UUIDs)
        # uuid_to_table : Dict[str, str] (UUID to table name)

        # make lazies
        logger.debug("Identifying classes for %d lazy proxies",
                     len(lazy_uuids))
        lazy_uuid_rows = self.backend.load_uuids_table(lazy_uuids)
        lazies = tools.group_by_function(lazy_uuid_rows,
                                         self.backend.uuid_row_to_table_name)
        new_uuids = self.proxy_factory.make_all_lazies(lazies)

        # get order and deserialize
        uuid_to_table_row = {r.uuid: r for r in to_load}
        ordered_uuids = get_reload_order(to_load, dependencies)
        new_uuids = self.deserialize_uuids(ordered_uuids, uuid_to_table,
                                           uuid_to_table_row, new_uuids)

        self.cache.update(new_uuids)
        results.update(new_uuids)
        new_results = [results[uuid] for uuid in input_uuids]

        # handle special case of storable functions
        for result in new_results:
            if isinstance(result, StorableFunction):
                self._sf_handler.register_function(result)

        return new_results


    def deserialize_uuids(self, ordered_uuids, uuid_to_table,
                          uuid_to_table_row, new_uuids=None):
        # TODO: remove this, replace with SerializationSchema
        logger.debug("Reconstructing from %d objects", len(ordered_uuids))
        new_uuids = tools.none_to_default(new_uuids, {})
        for uuid in ordered_uuids:
            if uuid not in self.cache and uuid not in new_uuids:
                # is_in = [k for (k, v) in dependencies.items() if v==uuid]
                table = uuid_to_table[uuid]
                table_row = uuid_to_table_row[uuid]
                table_dict = {attr: getattr(table_row, attr)
                              for (attr, type_name) in self.schema[table]}
                deserialize = self.class_info[table].deserializer
                obj = deserialize(uuid, table_dict, [new_uuids, self.cache])
                new_uuids[uuid] = obj
        return new_uuids

    def sync(self):
        pass

    def sync_all(self):
        pass

    def _cache_simulation_objects(self):
        # load up all the simulation objects
        try:
            backend_iter = itertools.chain(
                self.backend.table_iterator('simulation_objects'),
                self.backend.table_iterator('storable_functions')
            )
            sim_obj_uuids = [row.uuid for row in backend_iter]
        except KeyError:
            # TODO: this should probably be a custom error; don't rely on
            # the error type this backend raises
            # this happens if no simulation objects are given in the
            # schema... there's not technically required
            objs = []
        else:
            objs = self.load(sim_obj_uuids)
        return {get_uuid(obj): obj for obj in objs}

    def _reset_fixed_cache(self):
        self.cache.fixed_cache = {}
        self.cache.fixed_cache.update(self._simulation_objects)
        self.cache.fixed_cache.update(self._sf_handler.canonical_functions)

    def _update_pseudo_tables(self, simulation_objects):
        # TODO: replace the pseudo_tables code here with a class
        for uuid, obj in simulation_objects.items():
            my_cls = None
            for (key, cls) in self.simulation_classes.items():
                if isinstance(obj, cls):
                    self._pseudo_tables[key].append(obj)
                    my_cls = cls
                    # self._pseudo_tables[key][uuid] = obj
                    # if obj.is_named:
                        # self._pseudo_tables[key][obj.name] = obj
            if my_cls is None:
                self._pseudo_tables['misc_simulation'].append(obj)

    def summary(self, detailed=False):
        """Return a string summary of this storage file.

        Parameters
        ----------
        detailed : bool
            whether to return the detailed description, where the simulation
            objects are divided into their various pseudotables
        """
        out_str = "File: " + self.backend.filename + "\n"
        # TODO: add size to the first line
        out_str += "Includes tables:\n"
        storage_tables = dict(self._storage_tables)  # make a copy
        if detailed:
            storage_tables.pop('simulation_objects')
            storage_tables.update(self._pseudo_tables)

        for (name, table) in storage_tables.items():
            out_str += "* " + name + ": " + str(len(table)) + " items\n"
        return out_str

    def __getattr__(self, attr):
        # override getattr to create iterators over the tables (stores)
        if attr in self._storage_tables:
            return self._storage_tables[attr]
        elif attr in self._pseudo_tables:
            return self._pseudo_tables[attr]
        else:
            raise AttributeError("'{}' object has no attribute '{}'"\
                                 .format(self.__class__.__name__, attr))


class MixedCache(abc.MutableMapping):
    """Combine a frozen cache and a mutable cache"""
    # TODO: benchmark with single dict instead; might be just as fast!
    def __init__(self, fixed_cache=None):
        self.fixed_cache = tools.none_to_default(fixed_cache, default={})
        self.cache = {}

    def clear(self):
        self.cache = {}

    def delete_items(self, list_of_items, error_if_missing=False):
        for item in list_of_items:
            if item in self:
                del self[item]
            elif error_if_missing:
                raise KeyError()  # TODO: message and check error type

    def reproxy(self, schema):
        # TODO: idea: turn loaded objects back into None for any proxy
        # objects in cache. This frees those things up for garbage
        # collection.
        pass

    def __getitem__(self, key):
        if key in self.fixed_cache:
            value = self.fixed_cache[key]
        else:
            value = self.cache[key]
        return value

    def __setitem__(self, key, value):
        self.cache[key] = value

    def __delitem__(self, key):
        try:
            del self.cache[key]
        except KeyError as err:
            if key in self.fixed_cache:
                raise TypeError("Can't delete from fixed cache")
            else:
                raise err

    def __len__(self):
        return len(self.fixed_cache) + len(self.cache)

    def __iter__(self):
        return itertools.chain(self.fixed_cache, self.cache)

class StorageTable(abc.Sequence):
    # NOTE: currently you still need to be able to hold the whole table in
    # memory ... at least, with the SQL backend.
    def __init__(self, storage, table):
        self.storage = storage
        self.table = table
        self.clear_cache_block_freq = 100
        self.iter_block_size = 100  # TODO: base it on the size of an object

    def __iter__(self):
        # TODO: ensure that this gives us things in idx order
        backend_iter = self.storage.backend.table_iterator(self.table)
        enum_iter = enumerate(tools.grouper(backend_iter,
                                            self.iter_block_size))
        for block_num, block in enum_iter:
            row_uuids = [row.uuid for row in block]
            loaded = self.storage.load(row_uuids)
            if block_num % self.clear_cache_block_freq == 0 and block_num != 0:
                self.storage.cache.clear()
            for obj in loaded:
                yield obj

    def __getitem__(self, item):
        len_self = len(self)
        if not (-len_self <= item < len_self):
            raise IndexError("table index out of range")

        if item < 0:
            item += len(self)

        row = self.storage.backend.table_get_item(self.table, item)
        return self.storage.load([row.uuid])[0]

    def __len__(self):
        return self.storage.backend.table_len(self.table)

    def cache_all(self):
        old_blocksize = self.iter_block_size
        self.iter_block_size = len(self)
        _ = list(iter(self))

    def save(self, obj):
        # this is to match with the netcdfplus API
        self.storage.save(obj)

    # TODO: subclass for MCSteps with additional method .ordered, returning
    # things in the order of the mccycle number -- also, manage special
    # caching

class PseudoTable(abc.MutableSequence):
    # TODO: use this in the main code
    # NOTE: This will require that the storage class append to it
    """List of objects that can be retrieved by index or name.

    PseudoTables are used to group simulation objects together.
    """
    def __init__(self, sequence=None):
        self._sequence = []
        self._uuid_to_obj = {}
        self._name_to_uuid = {}
        sequence = tools.none_to_default(sequence, [])
        for item in sequence:
            self.append(item)

    @staticmethod
    def _get_uuid_and_name(obj):
        uuid = get_uuid(obj)
        try:
            name = None if not obj.is_named else obj.name
        except AttributeError:
            # occurs if simulation object is not a StorableNamedObject
            # (relevant for a few very old classes; should be fixed in 2.0)
            name = None
        return uuid, name

    def get_by_uuid(self, uuid):
        return self._uuid_to_obj[uuid]

    # NOTE: .index() can get confusing because you can have two equal
    # volumes (same CV, same range) with one named and the other not named.

    def __getitem__(self, item):
        try:
            ret_val = self._sequence[item]
        except TypeError:
            uuid = self._name_to_uuid[item]
            ret_val = self.get_by_uuid(uuid)
        return ret_val

    def __setitem__(self, key, value):
        # TODO: should this be allowed? or make it not really mutable, only
        # appendable?
        del self[key]
        self.insert(key, value)

    def __delitem__(self, key):
        # TODO: should this be allowed? or make it not really mutable, only
        # appendable?
        item = self[key]
        uuid, name = self._get_uuid_and_name(item)
        del self._sequence[self._sequence.index(item)]
        del self._uuid_to_obj[uuid]
        del self._name_to_uuid[name]

    def __len__(self):
        return len(self._uuid_to_obj)

    def __iter__(self):
        return iter(self._uuid_to_obj.values())

    def insert(self, where, item):
        # TODO: should this be allowed? or make it not really mutable, only
        # appendable?
        uuid, name = self._get_uuid_and_name(item)
        self._sequence.insert(where, item)
        self._uuid_to_obj[uuid] = item
        if name is not None:
            self._name_to_uuid[name] = uuid
