import os
import collections
import sqlalchemy as sql

# dict to convert from OPS string type descriptors to SQL types
sql_type = {
    'uuid': sql.String,
    'str': sql.String,
    'json': sql.String,
    'int': sql.Integer  #TODO add more
}

class SQLStorageBackend(object):
    """Generic storage backend for SQL.

    Uses SQLAlchemy; could easily duck-type an object that implements the
    necessary methods for other backends.
    """
    def __init__(self, filename, mode='r', template=None, fallback=None,
                 backend=None):
        self.template = template
        self.filename = filename
        self.fallback = fallback
        self.mode = mode
        self._metadata = sql.MetaData()
        if backend is None:
            backend = 'sqlite'

        if self.mode == "w" and os.path.exists(filename):
            # delete existing file; write after
            os.remove(filename)

        # we prevent writes by disallowing write method in read mode;
        # for everything else; just connect to the database
        connection_uri = self.filename_from_backend(filename, backend)
        self.engine = sql.create_engine(connection_uri)
        self.schema = {}

    @property
    def metadata(self):
        return self._metadata

    @property
    def class_to_table(self):
        # dict that links class to table it is stored in
        pass

    @staticmethod
    def filename_from_backend(filename, backend):
        # take backends like "sqlite", etc and return proper file connection
        # URI; would be essentially no-op for regular file opening as with
        # .nc
        uri_root = {
            'sqlite': "sqlite:///{filename}",
        }[backend]
        return uri_root.format(filename=filename)

    @staticmethod
    def _extract_metadata(sql_schema_metadata, table, column):
        if sql_schema_metadata:
            try:
                table_metadata = sql_schema_metadata[table]
            except KeyError:
                return {}
            try:
                col_metadata = table_metadata[column]
            except KeyError:
                return {}
            return col_metadata
        else:
            return {}

    ### FROM HERE IS THE GENERIC PUBLIC API
    def register_schema(self, schema, sql_schema_metadata=None):
        """Register (part of) a schema (create necessary tables in DB)

        Raises
        ------
        TypeError
            if the schema provided has names in the existing schema (may be
            trying to modify existing schema)

        """
        for table in schema:
            columns = [
                sql.Column(
                    col, sql_type[type_name],
                    **self._extract_metadata(sql_schema_metadata,
                                             table, col)
                )
                for (col, type_name) in schema[table]
            ]
            try:
                table = sql.Table(table, self.metadata, *columns)
            except sql.exc.InvalidRequestError:
                raise TypeError("Schema registration problem. Your schema "
                                "may already have tables of the same names.")

        self.metadata.create_all(self.engine)
        self.schema.update(schema)


    def add_to_table(self, objects):
        """Add a list of objects of a given class"""
        table = self.class_to_table[obj_class]
        # this will insert objects into the table
        pass

    def load_table_data(self, uuids):
        # this pulls out a table the information for the relevant UUIDs
        pass


class GeneralStorage(object):
    def __init__(self, filename, mode='r', template=None, fallback=None,
                 backend=None):
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


    def _cache_simulation_objects(self):
        # load up all the simulation objects
        pass

    def _create_virtual_stores(self, store_categories):
        # create virtual stores for simulation objects (e.g., .volume, etc)
        pass

    def save(self, obj):
        # prepare a single object for storage
        pass

    def save_list(self, list_of_objs):
        self.db.add_to_table(list_of_objs)

ops_schema = {}
ops_schema_sql_metadata = {}


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

