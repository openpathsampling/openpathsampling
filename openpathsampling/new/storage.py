import os
import collections
import openpathsampling as paths

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

class LazyLoader(object):
    def __init__(self, uuid, storage):
        self.uuid = self
        self.storage = storage
        self._loaded_object = None

    def __getattr__(self, attr):
        if self._loaded_object is None:
            self._loaded_object = storage.load(self.uuid, lazy=False)
        return getattr(self._loaded_object, attr)

    def repr(self):
        return "<LazyLoader for " + str(self.uuid) + ">"


class GeneralStorage(object):
    def __init__(self, filename, mode='r', template=None, fallback=None,
                 backend=None):

        self._init_new()
        self.simulation_objects = self._cache_simulation_objects()

    def _init_new(self):
        """Initial empty version of various dictionaries"""
        self.schema = {}
        self.serialization = {}
        self.class_to_table = {}

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

    def register_schema(self, schema, backend_metadata, class_to_table,
                        serialization):
        # check validity
        self.backend.register_schema(schema, backend_metadata)
        self.schema.update(schema)
        self.class_to_table.update(class_to_table)
        self.serialization.update(serialization)
        pass

    def table_for_class(self, class_):
        return self.class_to_table[class_]

    def _create_virtual_stores(self, store_categories):
        # create virtual stores for simulation objects (e.g., .volume, etc)
        pass

    def load(self, uuid, lazy=True):
        pass

    def save(self, obj):
        # prepare a single object for storage
        # check if obj is in DB
        pass

    def save_list(self, list_of_objs):
        # get the UUIDs of all objs
        # figure out which objs are already in the DB
        # organize by type
        # convert object to appropriate dict
        # gather UUIDs to construct
        pass

ops_schema = {
    'samples': {},
    'sample_sets': {},
    'trajectories': {},
    'move_changes': {},
    'steps': {}
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

class StorageCache(object):
    def __init__(self, n_subcaches):
        pass

class StorageList(object):
    def __init__(self):
        # set self.cache
        pass

    def __iter__(self):
        # iter fills the cache
        pass

    def __getitem__(self):
        # getitem builds the complete object
        # is is possible that using local lists of UUIDs to get might make
        # this just as fast? (stopping at trajectory level; no snapshots)
        pass

class SampleStore(StorageList):
    def cache_table(self, start_idx, end_idx):
        pass

    def cache_table_to_husks(self, cache_table):
        pass

    def fill_husks(self, husks):
        pass

