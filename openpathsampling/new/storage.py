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

Column names:

* ``uuid``
* ``idx``
"""

# these two tables are required in *all* schema
universal_schema = {
    'uuid': [('uuid', 'uuid'), ('table', 'int'), ('row', 'int')],
    'tables': [('name', 'str'), ('idx', 'int')]
}

class GeneralStorage(object):
    def __init__(self, filename, mode='r', template=None, fallback=None,
                 backend=None):
        self.simulation_objects = self._cache_simulation_objects()
        pass

    @classmethod
    def from_backend(cls, backend):
        obj = cls.__new__()
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

    def _create_virtual_stores(self, store_categories):
        # create virtual stores for simulation objects (e.g., .volume, etc)
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

ops_schema = {}
ops_schema_sql_metadata = {}
ops_class_to_table = {
    # this is only for data objects
    paths.Sample: 'samples',
    paths.SampleSet: 'sample_sets',
    paths.Trajectory: 'trajectories',
    paths.MoveChange: 'move_changes'
}

ops_class_to_serialization = {
    paths.Sample: (paths.Sample.to_dict, paths.Sample.from_dict),
    # this allows us to override the to_dict behavior for new storage
}

class OPSStorage(GeneralStorage):
    pass

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

