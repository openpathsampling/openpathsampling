import os
import collections
import itertools
import openpathsampling as paths
from serialization_helpers import get_uuid, get_all_uuids
from serialization_helpers import to_json_obj as serialize_sim
from serialization_helpers import from_json_obj as deserialize_sim
from serialization_helpers import get_all_uuids_loading
from serialization_helpers import dependency_dag, dag_reload_order
# TODO: both of these are from
from serialization_helpers import to_dict_with_uuids as serialize_data
from serialization_helpers import deserialize as deserialize_data
from class_info import ClassInfo, ClassInfoContainer
import tools
from openpathsampling.netcdfplus import StorableObject
from serialization import Serialization, DefaultDeserializer

import logging
logger = logging.getLogger(__name__)

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

# these two tables are required in *all* schema
universal_schema = {
    'uuid': [('uuid', 'uuid'), ('table', 'int'), ('row', 'int')],
    'tables': [('name', 'str'), ('idx', 'int'), ('module', 'str'),
               ('class_name', 'str')]
}

class GeneralStorage(object):
    def __init__(self, backend, class_info, schema=None, fallbacks=None):
        self.backend = backend
        self.schema = schema
        self.class_info = class_info
        self.mode = self.backend.mode
        # TODO: implement fallbacks
        self.fallbacks = tools.none_to_default(fallbacks, [])
        self.simulation_objects = self._cache_simulation_objects()
        self.cache = MixedCache(self.simulation_objects)
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
        logger.info("Missing info from {} dynamically-registered tables"\
                    .format(n_missing))
        classes = [table_to_class[tbl] for tbl in missing_info_tables]
        self.register_from_tables(missing_info_tables, classes)
        missing_info_tables = [tbl for tbl in self.schema
                               if tbl not in self.class_info.tables]
        logger.info("Successfully registered {} missing tables"\
                    .format(n_missing - len(missing_info_tables)))

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
        for info in class_info_list:
            info.set_defaults(schema)
            self.class_info.add_class_info(info)

        table_to_class = {table: self.class_info[table].cls
                          for table in schema
                          if table not in ['uuid', 'tables']}
        self.backend.register_schema(schema, table_to_class,
                                     backend_metadata)
        self.schema.update(schema)
            # here's where we add the class_info to the backend
        self.serialization.register_serialization(schema, self.class_info)


    def register_from_instance(self, lookup, obj):
        raise NotImplementedError("No way to register from an instance")

    def register_missing_tables_for_objects(self, uuid_obj_dict):
        # mistting items are handled by the special_lookup
        lookup_examples = {}
        for obj in uuid_obj_dict.values():
            lookup = self.class_info.lookup_key(obj)
            if lookup not in lookup_examples:
                self.register_from_instance(lookup, obj)

    def save(self, obj):
        # check if obj is in DB (maybe this can be removed?)
        exists = self.backend.load_uuids_table(uuids=[get_uuid(obj)],
                                               ignore_missing=True)
        if exists:
            return
        # find all UUIDs we need to save with this object
        # TODO: (perf) is this faster if we stop traversal on cached UUID?
        uuids = get_all_uuids(obj)
        # remove any UUIDs that have already been saved
        exists = self.backend.load_uuids_table(uuids=list(uuids.keys()),
                                               ignore_missing=True)
        for existing in exists:
            del uuids[existing.uuid]
        # group by table, then save appropriately
        # by_table; convert a dict of {uuid: obj} to {table: {uuid: obj}}
        get_table_name = lambda uuid, obj_: \
                self.class_info[obj_].table

        by_table = tools.dict_group_by(uuids, key_extract=get_table_name)

        # check default table for things to register; register them
        # TODO: convert to while?
        if '__missing__' in by_table:
            # __missing__ is a special result returned by the
            # ClassInfoContainer if this is object is expected to have a
            # table, but the table doesn't exist (e.g., for dynamically
            # added tables)
            missing = by_table.pop('__missing__')
            print missing
            logger.info("Attempting to register for {} missing objects".\
                        format(len(missing)))
            self.register_missing_tables_for_objects(missing)
            missing_by_table = tools.dict_group_by(missing, get_table_name)
            logger.info("Registered {} new tables: {}".\
                        format(len(missing_by_table),
                               str(list(missing_by_table.keys()))))
            by_table.update(missing_by_table)

        # this is the actual serialization
        for table in by_table:
            storables_list = [self.serialization.serialize[table](o)
                              for o in by_table[table].values()]
            logger.info("Storing {} objects to table {}".\
                        format(len(storables_list), table))
            self.backend.add_to_table(table, storables_list)
            logger.info("Storing complete")

    def load(self, input_uuids):
        # loading happens in 4 parts:
        # 1. Get UUIDs that need to be loaded
        # 2. Build the DAG to determine loading order
        # 3. Make all lazy objects
        # 4. Reserialize remaining in DAG order
        results = {uuid: self.cache[uuid] for uuid in input_uuids
                   if uuid in self.cache}
        uuid_list = [uuid for uuid in input_uuids if uuid not in self.cache]
        to_load, lazies, dependencies, uuid_to_table = \
                get_all_uuids_loading(uuid_list=uuid_list,
                                      backend=self.backend,
                                      schema=self.schema,
                                      existing_uuids=self.cache)

        # build the dag
        dag = dependency_dag(dependencies)
        # lazies can't have dependencies
        lazies = tools.group_by_function(lazies,
                                         lambda r: uuid_to_table[r.uuid])
        new_uuids = self.serialization.make_all_lazies(lazies)


        # deserialize in order
        uuid_to_table_row = {r.uuid: r for r in to_load}
        for uuid in dag_reload_order(dag):
            table = uuid_to_table[uuid]
            table_row = uuid_to_table_row[uuid]
            table_dict = {attr: getattr(table_row, attr)
                          for (attr, type_name) in self.schema[table]}
            # TODO: improve this
            deserialize = self.class_info[table].deserializer
            obj = deserialize(uuid, table_dict, [new_uuids, self.cache])
            new_uuids[uuid] = obj

        self.cache.update(new_uuids)
        results.update(new_uuids)
        return [results[uuid] for uuid in input_uuids]

    def _cache_simulation_objects(self):
        # load up all the simulation objects
        return {}

    def _create_virtual_stores(self, store_categories):
        # create virtual stores for simulation objects (e.g., .volume, etc)
        pass

    # def __getattr__(self, attr):
        # override getattr to create iterators over the tables (stores)
        # if attr in schema:
            # return TableIterator(self, attr)


class MixedCache(collections.MutableMapping):
    """Combine a frozen cache and a mutable cache"""
    # TODO: benchmark with single dict instead; might be just as fast!
    def __init__(self, fixed_cache=None):
        self.fixed_cache = tools.none_to_default(fixed_cache, default={})
        self.cache = {}

    def __getitem__(self, key):
        if key in self.fixed_cache:
            return self.fixed_cache[key]
        else:
            return self.cache[key]

    def __setitem__(self, key, value):
        self.cache[key] = value

    def __delitem__(self, key, value):
        try:
            del self.cache[key]
        except KeyError as e:
            if key in self.fixed_cache:
                raise TypeError("Can't delete from fixed cache")
            else:
                raise e

    def __len__(self):
        return len(self.fixed_cache) + len(self.cache)

    def __iter__(self):
        return itertools.chain(self.fixed_cache, self.cache)


# TODO: collections.Sequence?
class StorageTable(object):
    def __init__(self, storage, table):
        self.storage = storage
        self.table = table

    def __iter__(self):
        backend_iterator = self.storage.backend.table_iterator(self.table)
        for row in backend_iterator:
            yield self.storage.load([row.uuid])

    def _get_single(self, value):
        pass

    def __getitem__(self):
        # getitem builds the complete object
        # is is possible that using local lists of UUIDs to get might make
        # this just as fast? (stopping at trajectory level; no snapshots)
        pass

    def store_order(self):
        # return in order of the idx
        pass

    # TODO: subclass for MCSteps with additional method .ordered, returning
    # things in the order of the mccycle number