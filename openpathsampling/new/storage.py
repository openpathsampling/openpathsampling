import os
import collections
import itertools
import openpathsampling as paths
from serialization import get_uuid, get_all_uuids
from serialization import to_json as serialize_sim
from serialization import from_json as deserialize_sim
# TODO: both of these are from
from serialization import to_dict_with_uuids as serialize_data
from serialization import deserialize as deserialize_data
import tools
from openpathsampling.netcdfplus import StorableObject

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

ClassInfo = collections.namedtuple(
    'ClassInfo',
    ['table', 'cls', 'serializer', 'deserializer']
)

#TODO: should this become a collections.Sequence?
class ClassInfoContainer(object):
    def __init__(self, default_info, class_info_list=None):
        class_info_list = tools.none_to_default(class_info_list, [])
        self.class_to_info = {}
        self.table_to_info = {}
        self.class_to_table = {}
        self.table_to_class = {}
        self.class_info_list = []
        self.default_info = default_info
        for info in class_info_list:
            self.add_class_info(info)

    def add_class_info(self, info_node):
        # check that we're not in any existing
        self.class_info_list.append(info_node)
        self.class_to_info.update({info_node.cls: info_node})
        self.table_to_info.update({info_node.table: info_node})
        self.class_to_table.update({info_node.cls: info_node.table})
        self.table_to_class.update({info_node.table: info_node.cls})

    def __getitem__(self, item):
        if tools.is_string(item):
            return self.table_to_info[item]
        else:
            try:
                return self.class_to_info[item]
            except KeyError as e:
                if issubclass(item, self.default_info.cls):
                    return self.default_info
                else:
                    raise e

    def __repr__(self):  # pragma: no cover
        return ("ClassInfoContainer(default_info=" + repr(self.default_info)
                + ", class_info_list=" + repr(self.class_info_list) + ")")



def make_lazy_class(cls_):
    # this is to mix-in inheritence
    class LazyClass(LazyLoader, cls_):
        pass
    return LazyClass

class LazyLoader(object):
    def __init__(self, uuid, class_, storage):
        self.__uuid__ = uuid
        self.storage = storage
        self.class_= class_
        self._loaded_object = None

    def load(self):
        if self._loaded_object is None:
            self._loaded_object = self.storage.load(self.__uuid__, lazy=False)
        return self._loaded_object

    def __getattr__(self, attr):
        return getattr(self.load(), attr)

    def __iter__(self):
        loaded = self.load()
        if not hasattr(loaded, '__iter__'):
            raise TypeError() # TODO: message
        # TODO: load all objects in the list
        return loaded.__iter__

    # TODO: how to handle isinstance? each lazy-loading class needs a sub
    def repr(self):
        return ("<LazyLoader for " + str(self.class_.__name__) + " UUID "
                + str(self.__uuid__) + ">")


class MixedCache(collections.MutableMapping):
    """Combine a frozen cache and a mutable cache"""
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
    def __init__(self, backend, schema, class_info, fallbacks=None):
        self.backend = backend
        self.schema = schema
        self.class_info = class_info
        self._lazy_classes = {}
        self.simulation_objects = self._cache_simulation_objects()
        self.cache = MixedCache(self.simulation_objects)
        self.register_schema(self.schema, class_info_list=[])

    def _cache_simulation_objects(self):
        # load up all the simulation objects
        return {}

    def make_lazy(self, table, uuid):
        class_ = self.table_to_class[table]
        if table not in self._lazy_classes:
            self._lazy_classes[table] = make_lazy_class(class_)
        return self._lazy_classes[table](uuid=uuid,
                                         cls=class_,
                                         storage=self)

    def register_schema(self, schema, class_info_list,
                        backend_metadata=None):
        # check validity
        self.backend.register_schema(schema, backend_metadata)
        self.schema.update(schema)
        for info in class_info_list:
            self.class_info.add_class_info(info)


    def table_for_class(self, class_):
        return self.class_to_table[class_]

    def _create_virtual_stores(self, store_categories):
        # create virtual stores for simulation objects (e.g., .volume, etc)
        pass

    def load(self, uuid, lazy=None):
        # get UUIDs and tables associated
        # if lazy, return the lazy object
        # if table has custom loader, use that
        pass

    def serialize(self, obj):
        class_info = self.class_info[obj.__class__]
        return class_info.serializer(obj)

    def save(self, obj):
        # check if obj is in DB
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
                self.class_info[obj_.__class__].table

        by_table = tools.dict_group_by(uuids, key_extract=get_table_name)

        for table in by_table:
            storables_list = [self.serialize(o)
                              for o in by_table[table].values()]
            self.backend.add_to_table(table, storables_list)


    def __getattr__(self, attr):
        # override getattr to create iterators over the tables
        pass

ops_schema = {
    'samples': [('trajectory', 'lazy'), ('ensemble', 'uuid'),
                ('replica', 'int'),
                # in my opinion, the next 3 should be removed
                ('parent', 'lazy'), ('bias', 'float'),
                ('mover', 'uuid')],
    'sample_sets': [('samples', 'list_uuid'), ('movepath', 'lazy')],
    'trajectories': [('snapshots', 'list_uuid')],
    'move_changes': [('mover', 'uuid'), ('details', 'lazy'), ('cls', 'str'),
                     ('subchanges', 'list_uuid'), ('samples', 'list_uuid'),
                     ('input_samples', 'list_uuid')],
    'steps': [('change', 'uuid'), ('active', 'uuid'), ('previous', 'uuid'),
              ('simulation', 'uuid'), ('mccycle', 'int')],
    'details': [('json', 'json')],
    'simulation_objects': [('json', 'json'), ('class_idx', 'int')]
}

ops_schema_sql_metadata = {}

ops_class_info = ClassInfoContainer(
    default_info=ClassInfo('simulation_objects', cls=StorableObject,
                           serializer=serialize_sim,
                           deserializer=deserialize_sim),
    class_info_list=[
        ClassInfo(table='samples', cls=paths.Sample,
                  serializer=serialize_data,
                  deserializer=deserialize_data),
        ClassInfo(table='sample_sets', cls=paths.SampleSet,
                  serializer=serialize_data,
                  deserializer=deserialize_data),
        ClassInfo(table='trajectories', cls=paths.Trajectory,
                  serializer=serialize_data,
                  deserializer=deserialize_data),
        ClassInfo(table='move_changes', cls=paths.MoveChange,
                  serializer=deserialize_data,
                  deserializer=deserialize_data),  #TODO: may need custoom
        ClassInfo(table='steps', cls=paths.MCStep,
                  serializer=serialize_data,
                  deserializer=deserialize_data)
    ]
)

class TableIterator(object):
    def __init__(self, storage):
        self.cache = storage.simulation_objects.copy()

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
