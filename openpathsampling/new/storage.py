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
            self._loaded_object = storage.load(self.__uuid__, lazy=False)
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
        return self.lazy_classes[table](uuid=uuid,
                                        cls=self.table_to_class[table],
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
              ('simulation', 'uuid'), ('mccycle', 'int')]
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
