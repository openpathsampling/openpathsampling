import os
import collections
import itertools
import openpathsampling as paths
from serialization_helpers import get_uuid, get_all_uuids
from serialization_helpers import to_json_obj as serialize_sim
from serialization_helpers import from_json_obj as deserialize_sim
# TODO: both of these are from
from serialization_helpers import to_dict_with_uuids as serialize_data
from serialization_helpers import deserialize as deserialize_data
import tools
from openpathsampling.netcdfplus import StorableObject
from serialization import Serialization

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
    'tables': [('name', 'str'), ('idx', 'int')]
}

class ClassInfo(object):
    # I so wanted this to be a namedtuple, but it needs more attention
    """
    Parameters
    ----------
    table : str
        the table name for this object type
    cls : class
        the class for this object type (used in the default deserializer)
    serializer : callable
        serializer for this object type (None gives the default serializer)
    deserializer : callable
        deserializer for this object type (None gives the default
        deserializer)
    lookup_result : any
        the result when ClassInfoContainer looks up objects in this table
    """
    def __init__(self, table, cls, serializer=None, deserializer=None,
                 lookup_result=None):
        self.table = table
        self.cls = cls
        self.serializer = serializer
        self.deserializer = deserializer
        self.lookup_result = lookup_result

class ClassInfoContainer(object):
    """Connect table to serialization information.

    This is the primary location for information connecting application
    objects (i.e., OPS) to the serialization mechanisms. In particular, the
    __getattr__ in this can take either an instance or a table name, and in
    either case links it to the ClassInfo. This lets us link an object
    instance to the name of the table where it should be stored, and to
    link the table to the serialization methods.

    Notes
    -----
    Linking instances to tables can require special treatment. The default
    behavior is to use the class of the object to identify the table, and
    this works for most cases. However, sometimes a class can map to more
    than one table, and therefore is not a unique key. For example, a single
    class might be used to represent data with different dimensions, and
    therefore require different tables (e.g, coordinates for different
    systems). In such cases, the ClassInfoContainer needs to be subclassed
    with specialized information.
    """
    def __init__(self, default_info, class_info_list=None):
        class_info_list = tools.none_to_default(class_info_list, [])
        self.lookup_to_info = {}
        self.table_to_info = {}
        self.class_info_list = []
        self.default_info = default_info
        self.missing_table = ClassInfo(table="__missing__", cls=None)
        self.add_class_info(default_info)
        for info in class_info_list:
            self.add_class_info(info)

    def add_class_info(self, info_node):
        # check that we're not in any existing
        self.class_info_list.append(info_node)
        # TODO: is it possible to generalize the special cases (use the type
        # of the expected return to identify the function to call -- won't
        # be completely general, but may simplify some)
        self.lookup_to_info.update({info_node.lookup_result: info_node})
        self.table_to_info.update({info_node.table: info_node})

    def lookup_key(self, item):
        if not self.is_special(item):
            lookup_key = item.__class__
        else:
            lookup_key = self.special_lookup_key(item)
        return lookup_key

    def is_special(self, item):
        return False

    def special_lookup_key(self, item):
        return NotImplementedError("No special types implemented")

    def get_special(self, item):
        lookup = self.special_lookup_key(item)
        try:
            return self.lookup_to_info[lookup]
        except KeyError:
            return self.missing_table

    def __getitem__(self, item):
        if tools.is_string(item):
            return self.table_to_info[item]
        elif self.is_special(item):
            return self.get_special(item)
        else:
            lookup = self.lookup_key(item)
            try:
                return self.lookup_to_info[lookup]
            except KeyError as e:
                if isinstance(item, self.default_info.cls):
                    return self.default_info
                else:
                    raise e

    def __repr__(self):  # pragma: no cover
        return ("ClassInfoContainer(default_info=" + repr(self.default_info)
                + ", class_info_list=" + repr(self.class_info_list) + ")")


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


class GeneralStorage(object):
    def __init__(self, backend, class_info, schema=None, fallbacks=None):
        self.backend = backend
        self.schema = schema
        self.class_info = class_info
        # TODO: implement fallbacks
        self.fallbacks = tools.none_to_default(fallbacks, [])
        self.simulation_objects = self._cache_simulation_objects()
        self.cache = MixedCache(self.simulation_objects)
        self.serialization = Serialization(self)
        if schema is None:
            schema = backend.schema
        self.register_schema(self.schema, class_info_list=[])

    def close(self):
        # TODO: should sync on close
        self.backend.close()
        for fallback in self.fallbacks:
            fallback.close()

    def register_schema(self, schema, class_info_list,
                        backend_metadata=None):
        # check validity
        self.backend.register_schema(schema, backend_metadata)
        self.schema.update(schema)
        for info in class_info_list:
            self.class_info.add_class_info(info)
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
            logger.info("Storing {} UUIDs to table {}".\
                        format(len(storables_list), table))
            self.backend.add_to_table(table, storables_list)
            logger.info("Storing complete")

    def load(self, uuid, lazy=None):
        # get UUIDs and tables associated
        # if lazy, return the lazy object
        # if table has custom loader, use that
        pass

    def _cache_simulation_objects(self):
        # load up all the simulation objects
        return {}

    def _create_virtual_stores(self, store_categories):
        # create virtual stores for simulation objects (e.g., .volume, etc)
        pass

    #def __getattr__(self, attr):
        # override getattr to create iterators over the tables (stores)
    #    pass


class TableIterator(object):
    def __init__(self, storage, table):
        self.storage = storage

    def __iter__(self):
        # iter manages the cache
        pass

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
