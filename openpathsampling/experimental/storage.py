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
import itertools

import tools
from serialization_helpers import get_uuid, get_all_uuids
from serialization_helpers import get_all_uuids_loading
from serialization_helpers import get_reload_order
from serialization import Serialization

logger = logging.getLogger(__name__)

# these two tables are required in *all* schema
universal_schema = {
    'uuid': [('uuid', 'uuid'), ('table', 'int'), ('row', 'int')],
    'tables': [('name', 'str'), ('idx', 'int'), ('module', 'str'),
               ('class_name', 'str')]
}

class GeneralStorage(object):
    def __init__(self, backend, class_info, schema=None,
                 simulation_classes=None, fallbacks=None):
        self.backend = backend
        self.schema = schema
        self.class_info = class_info
        self.mode = self.backend.mode
        # TODO: implement fallbacks
        self.fallbacks = tools.none_to_default(fallbacks, [])

        self.simulation_classes = tools.none_to_default(simulation_classes,
                                                        {})

        self._pseudo_tables = {table_name: dict()
                               for table_name in self.simulation_classes}

        self._storage_tables = {}  # stores .steps, .snapshots
        self._simulation_objects = self._cache_simulation_objects()
        self.cache = MixedCache(self._simulation_objects)
        self.serialization = Serialization(self)
        if self.schema is None:
            self.schema = backend.schema
        self.initialize_with_mode(self.mode)

    def initialize_with_mode(self, mode):
        if mode == 'r' or mode == 'a':
            table_to_class = self.backend.table_to_class
            self._load_missing_info_tables(table_to_class)

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

    def close(self):
        # TODO: should sync on close
        self.backend.close()
        for fallback in self.fallbacks:
            fallback.close()

    def register_schema(self, schema, class_info_list,
                        backend_metadata=None, read_mode=False):
        # check validity
        self.class_info.register_info(class_info_list, schema)
        # for info in class_info_list:
            # info.set_defaults(schema)
            # self.class_info.add_class_info(info)

        table_to_class = {table: self.class_info[table].cls
                          for table in schema
                          if table not in ['uuid', 'tables']}
        # here's where we add the class_info to the backend
        self.backend.register_schema(schema, table_to_class,
                                     backend_metadata)
        self.schema.update(schema)
        for table in self.schema:
            self._storage_tables[table] = StorageTable(self, table)
        self.serialization.register_serialization(schema, self.class_info)


    def register_from_instance(self, lookup, obj):
        raise NotImplementedError("No way to register from an instance")

    def register_missing_tables_for_objects(self, uuid_obj_dict):
        # mistting items are handled by the special_lookup
        lookup_examples = set([])
        for obj in uuid_obj_dict.values():
            lookup = self.class_info.lookup_key(obj)
            if lookup not in lookup_examples:
                self.register_from_instance(lookup, obj)
                lookup_examples |= {lookup}

    def uuids_in_storage(self, uuid_list):
        return self.backend.load_uuids_table(uuids=uuid_list,
                                             ignore_missing=True)

    def save(self, obj):
        # TODO: convert the whole .save process to something based on the
        # class_info.serialize method (enabling per-class approaches for
        # finding UUIDs, which will be a massive serialization speed-up
        # self.class_info.serialize(obj, storage=self)
        # check if obj is in DB (maybe this can be removed?)
        logger.debug("Starting save")
        exists = self.backend.load_uuids_table(uuids=[get_uuid(obj)],
                                               ignore_missing=True)
        if exists:
            return
        # find all UUIDs we need to save with this object
        logger.debug("Listing all objects to save")
        uuids = get_all_uuids(obj, known_uuids=self.cache,
                              class_info=self.class_info)
        logger.debug("Checking if objects already exist in database")
        # remove any UUIDs that have already been saved
        exists = self.backend.load_uuids_table(uuids=list(uuids.keys()),
                                               ignore_missing=True)
        for existing in exists:
            del uuids[existing.uuid]
        # group by table, then save appropriately
        # by_table; convert a dict of {uuid: obj} to {table: {uuid: obj}}
        get_table_name = lambda uuid, obj_: self.class_info[obj_].table

        by_table = tools.dict_group_by(uuids, key_extract=get_table_name)

        # check default table for things to register; register them
        # TODO: convert to while?
        if '__missing__' in by_table:
            # __missing__ is a special result returned by the
            # ClassInfoContainer if this is object is expected to have a
            # table, but the table doesn't exist (e.g., for dynamically
            # added tables)
            missing = by_table.pop('__missing__')
            logger.info("Registering tables for %d missing objects",
                        len(missing))
            self.register_missing_tables_for_objects(missing)
            missing_by_table = tools.dict_group_by(missing, get_table_name)
            logger.info("Registered %d new tables: %s",
                        len(missing_by_table),
                        str(list(missing_by_table.keys())))
            by_table.update(missing_by_table)

        # this is the actual serialization
        logger.debug("Filling %d tables: %s", len(by_table),
                     str(list(by_table.keys())))
        for table in by_table:
            logger.debug("Storing %d objects to table %s",
                         len(by_table[table]), table)
            serialize = self.class_info[table].serializer

            # DEBUG
            # if table == 'move_changes':
                # for o in by_table[table].values():
                    # print o
                    # serialize(o)

            storables_list = [serialize(o) for o in by_table[table].values()]
            self.backend.add_to_table(table, storables_list)
            # special handling for simulation objects
            if table == 'simulation_objects':
                self._update_pseudo_tables(by_table[table])
            logger.debug("Storing complete")

    def load(self, input_uuids, force=False):
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
                                      existing_uuids=self.cache)
        logger.debug("Loading %d objects; creating %d lazy proxies",
                     len(to_load), len(lazy_uuids))

        # make lazies
        logger.debug("Identifying classes for %d lazy proxies",
                     len(lazy_uuids))
        lazy_uuid_rows = self.backend.load_uuids_table(lazy_uuids)
        lazies = tools.group_by_function(lazy_uuid_rows,
                                         self.backend.uuid_row_to_table_name)
        # TODO: replace this with something not based on Serialization
        # object
        new_uuids = self.serialization.make_all_lazies(lazies)

        # get order and deserialize
        uuid_to_table_row = {r.uuid: r for r in to_load}
        ordered_uuids = get_reload_order(to_load, dependencies)
        new_uuids = self.deserialize_uuids(ordered_uuids, uuid_to_table,
                                           uuid_to_table_row, new_uuids)

        self.cache.update(new_uuids)
        results.update(new_uuids)
        return [results[uuid] for uuid in input_uuids]

    def deserialize_uuids(self, ordered_uuids, uuid_to_table,
                          uuid_to_table_row, new_uuids=None):
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
        # backend_iterator = self.backend.table_iterator('simulation_objects')
        # sim_obj_uuids = [row.uuid for row in backend_iterator]
        # objs = self.load(sim_obj_uuids)
        # load up all the simulation objects
        return {}

    def _update_pseudo_tables(self, simulation_objects):
        for uuid, obj in simulation_objects.items():
            for (key, cls) in self.simulation_classes.items():
                if isinstance(obj, cls):
                    self._pseudo_tables[key][uuid] = obj
                    if obj.is_named:
                        self._pseudo_tables[key][obj.name] = obj
                    continue

    def __getattr__(self, attr):
        # override getattr to create iterators over the tables (stores)
        if attr in self._storage_tables:
            return self._storage_tables[attr]
        elif attr in self._pseudo_tables:
            return self._pseudo_tables[attr]
        else:
            raise AttributeError("'{}' object has no attribute '{}'"\
                                 .format(self.__class__.__name__, attr))


class MixedCache(collections.MutableMapping):
    """Combine a frozen cache and a mutable cache"""
    # TODO: benchmark with single dict instead; might be just as fast!
    def __init__(self, fixed_cache=None):
        self.fixed_cache = tools.none_to_default(fixed_cache, default={})
        self.cache = {}

    def delete_items(self, list_of_items, error_if_missing=False):
        for item in list_of_items:
            if item in self:
                del self[item]
            elif error_if_missing:
                raise KeyError()  # TODO: message and check error type

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


class StorageTable(collections.Sequence):
    # NOTE: currently you still need to be able to hold the whole table in
    # memory ... at least, with the SQL backend.
    def __init__(self, storage, table, cache=None):
        self.storage = storage
        self.table = table
        self.clear_cache_frequency = 1
        self.iter_block_size = 1

    def __iter__(self):
        # TODO: ensure that this gives us things in idx order
        backend_iterator = self.storage.backend.table_iterator(self.table)
        # TODO: implement use of self.iter_block_size
        # TODO: implement use of self.clear_cache_frequency
        for row in backend_iterator:
            yield self.storage.load([row.uuid])[0]

    def __getitem__(self, item):
        backend_iterator = self.storage.backend.table_iterator(self.table)
        if item < 0:
            item += len(self)
        n_iter = 0
        row = backend_iterator.next()
        while row and n_iter < item:
            row = backend_iterator.next()
            n_iter += 1
        return self.storage.load([row.uuid])[0]

    def __len__(self):
        backend_iterator = self.storage.backend.table_iterator(self.table)
        return len(list(backend_iterator))

    def save(self, obj):
        # this is to match with the netcdfplus API
        self.storage.save(obj)

    # TODO: subclass for MCSteps with additional method .ordered, returning
    # things in the order of the mccycle number -- also, manage special
    # caching
