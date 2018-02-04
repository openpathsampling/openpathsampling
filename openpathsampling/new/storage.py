import os
import collections
import itertools
import openpathsampling as paths
from serialization import get_uuid
import tools

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
    ['table', 'class', 'serializer', 'deserializer']
)

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
    def __init__(self, filename, mode='r', template=None, fallback=None,
                 backend=None):
        self._init_new()
        self.simulation_objects = self._cache_simulation_objects()
        self.cache = MixedCache(self.simulation_objects)

    def _init_new(self):
        """Initial empty version of various dictionaries"""
        self.schema = {}
        self.serialization = {}
        self.class_to_table = {}
        self.table_to_class = {}
        self.lazy_tables = {}

    @classmethod
    def from_backend(cls, backend):
        obj = cls.__new__()
        obj._init_new()
        obj.filename = backend.filename
        obj.mode = backend.mode
        obj.template = backend._template
        obj.fallback = backend.fallback
        obj.backend = backend.backend
        obj.db = backend
        obj.simulation_objects = obj._cache_simulation_objects()

    def _cache_simulation_objects(self):
        # load up all the simulation objects
        pass

    def make_lazy(self, table, uuid):
        class_ = self.table_to_class[table]
        if table not in self.lazy_classes:
            self.lazy_classes[table] = make_lazy_class(class_)
        return self.lazy_classes[table](uuid=uuid,
                                        cls=class_,
                                        storage=self)

    def register_schema(self, schema, class_to_table, serialization,
                        lazy_tables=None, backend_metadata=None):
        # check validity
        self.backend.register_schema(schema, backend_metadata)
        self.schema.update(schema)
        self.class_to_table.update(class_to_table)
        self.table_to_class.update({
            table: cls_ for (cls_, table) in class_to_table.items()
        })
        self.serialization.update(serialization)
        if lazy_tables:
            if lazy_tables is True:
                lazy_tables = list(schema.keys())
            self.lazy_tables += lazy_tables

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
        class_info = self.get_class_info(obj.__class__)
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
                self.get_class_info(obj_.__class__).table

        by_table = dict_group_by(uuids, key_extract=get_table_name)

        for table in by_table:
            storables_list = [self.serialize(o) for o in by_table[table]]
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
    'details': [('json', 'json')]
    'simulation_objects': [('json', 'json'), ('class_idx', 'int')]
}
ops_schema_sql_metadata = {}
ops_class_to_table = {
    # this is only for data objects
    paths.Sample: 'samples',
    paths.SampleSet: 'sample_sets',
    paths.Trajectory: 'trajectories',
    paths.MoveChange: 'move_changes',
    paths.MCStep: 'steps'
}

ops_class_to_serialization = {
    paths.Sample: (paths.Sample.to_dict, paths.Sample.from_dict),
    # this allows us to override the to_dict behavior for new storage
}

class OPSStorage(GeneralStorage):
    def table_for_class(self, class_):
        try:
            table = ops_class_to_table[class_]
        except KeyError:
            if issubclass(class_, paths.netcdfplus.StorableNamedObject):
                table = 'simulation_object'
            else:
                table = None
        return table

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
