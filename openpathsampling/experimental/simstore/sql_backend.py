import os
import collections
from collections import abc
import sqlalchemy as sql
from .storage import universal_schema
from .tools import group_by, compare_sets, grouper
from . import tools
# import ujson as json  # ujson is no longer maintained
import json

from .backend import extract_backend_metadata
from .my_types import backend_registration_type
from .serialization_helpers import import_class


import logging
logger = logging.getLogger(__name__)

# dict to convert from OPS type descriptors to SQL types
sql_type = {
    'uuid': sql.String,
    'lazy': sql.String,
    'list_uuid': sql.String,
    'str': sql.String,
    'json': sql.String,
    'json_obj': sql.String,
    'int': sql.Integer,
    'float': sql.Float,
    'function': sql.String,
    'ndarray': sql.LargeBinary,
    'bool': sql.Boolean,
    #TODO add more
}

universal_sql_meta = {
    'uuid': {'uuid': {'primary_key': True}},
    'tables': {'name': {'primary_key': True}}
}

def make_columns(table_name, schema, sql_schema_metadata, backend_types):
    columns = []
    type_mapping = {k: v[0] for k, v in backend_types.items()}
    # TODO: use size_info for fixed-width columns
    size_info = {k: v[1] for k, v in backend_types.items()}
    if table_name not in universal_schema:
        columns.append(sql.Column('idx', sql.Integer,
                                  primary_key=True))
        columns.append(sql.Column('uuid', sql.String))
    for col, type_name in schema[table_name]:
        col_type = sql_type[type_mapping[type_name]]
        metadata = extract_backend_metadata(sql_schema_metadata,
                                            table_name, col)
        columns.append(sql.Column(col, col_type, **metadata))
    return columns


# TODO: test this, and then see if I can refactor existing code to use it
# NOTE: looks like the results have to be loaded into memory before this
# function ends (otherwise the DB connection closes and errors .. so make a
# list of them at the end, only use when we want the whole list?
def sql_db_execute(engine, statements, execmany_dict=None):
    """Convenience for interacting with the database."""
    if not isinstance(statements, list):
        listified = True
        statements = [statements]
    if execmany_dict is None:
        execmany_dict = {}
    with engine.connect() as conn:
        results = [conn.execute(statement, execmany_dict)
                   for statement in statements]
    if listified:
        results = results[0]
    return results

from openpathsampling.netcdfplus import StorableNamedObject
class SQLStorageBackend(StorableNamedObject):
    """Generic storage backend for SQL.

    Uses SQLAlchemy; could easily duck-type an object that implements the
    necessary methods for other backends.

    Parameters
    ----------
    filename : str
        name of the sqlite file or connection URI for a service-based
        database
    mode : 'r', 'w', or 'a'
        "file mode", as with files -- this is a little strange for
        databases, but keeps the interface consistent
    sql_dialect : str
        name of the SQL dialect to use; default is sqlite

    Additional keyword arguments are passed to sqlalchemy.create_engine. Of
    particular use is ``echo`` (bool) which echos SQL commands to stdout
    (useful for debugging).

    More info: https://docs.sqlalchemy.org/en/latest/core/engines.html
    """
    MAX_SQL_ITEMS = 900
    def __init__(self, filename, mode='r', sql_dialect='sqlite', **kwargs):
        super().__init__()
        self.filename = filename
        self.sql_dialect = sql_dialect
        self.mode = mode
        self.kwargs = kwargs
        self.debug = False

        # maps a specific type name, to generic type info, e.g.
        # 'ndarray.float32(1651,3)': 'ndarray'
        # keys here are in the schema, values are (sql_type, size)
        self.known_types = {k: (k, None) for k in sql_type}

        # override later if mode == 'r' or 'a'
        self.schema = {}
        self.table_to_number = {}
        self.number_to_table = {}

        # TODO: change this to an LRU cache
        self.known_uuids = collections.defaultdict(set)

        # we prevent writes by disallowing write method in read mode;
        # for everything else; just connect to the database
        self.engine = None
        self.connection_uri = None
        self._metadata = None
        if filename is not None:
            file_exists = os.path.exists(filename)
            if self.mode == "w" and file_exists:
                # delete existing file; write after
                os.remove(filename)

            if self.mode == 'a' and not file_exists:
                # act as if the mode is 'w'; note we change this back later
                self.mode = 'w'

            self.connection_uri = self.filename_from_dialect(
                filename,
                self.sql_dialect
            )
            engine = sql.create_engine(self.connection_uri, **self.kwargs)
            self._initialize_from_engine(engine)
            self.mode = mode  # in case we changed when checking existence

    def _initialize_from_engine(self, engine):
        self.engine = engine
        self._metadata = sql.MetaData(bind=self.engine)
        self._initialize_with_mode(self.mode)

    @property
    def identifier(self):
        if self.connection_uri == "sqlite:///:memory:":
            return "sqlite:///{}".format(id(self)), self.mode
        else:
            return self.connection_uri, self.mode

    def to_dict(self):
        return {
            'filename': self.filename,
            'sql_dialect': self.sql_dialect,
            'mode': 'a' if self.mode == 'w' else self.mode,
            'connection_uri': self.connection_uri,
            'kwargs': self.kwargs,
        }

    @classmethod
    def from_dict(cls, dct):
        init_dct = dct.copy()
        filename = init_dct.pop('filename')
        connection_uri = init_dct.pop('connection_uri')
        kwargs = init_dct.pop('kwargs')
        init_dct.update(kwargs)
        obj = cls(filename=None, **init_dct)
        obj.connection_uri = connection_uri
        obj.filename = filename
        engine = sql.create_engine(connection_uri, **obj.kwargs)
        obj._initialize_from_engine(engine)
        return obj

    def __reduce__(self):
        return (self.from_dict, (self.to_dict(),))

    def _initialize_with_mode(self, mode):
        """setup of tables; etc, as varies between w and r/a modes"""
        if mode == "w":
            # TODO: drop tables
            schema_table = sql.Table('schema', self.metadata,
                                     sql.Column('table', sql.String),
                                     sql.Column('schema', sql.String))
            metadata_table = sql.Table('metadata', self.metadata,
                                       sql.Column('key', sql.String),
                                       sql.Column('value', sql.String))

            self.metadata.create_all(self.engine)
            self.register_schema(universal_schema, universal_sql_meta)
            uuid_table = self.metadata.tables['uuid']
            index = sql.Index('uuids_index', *uuid_table.c, unique=True)
            index.create(self.engine)

        elif mode == "r" or mode == "a":
            self.metadata.reflect(self.engine)
            self.schema = self.database_schema()
            self.table_to_number, self.number_to_table = \
                    self.internal_tables_from_db()

    @classmethod
    def from_engine(cls, engine, connection_uri=None, **kwargs):
        """Constructor allowing user to specify the SQLAlchemy Engine.

        The default constructor doesn't allow all the options to the engine
        that SQLAlchemy can provide. Use this if you want more customization
        of the database engine.

        More info: https://docs.sqlalchemy.org/en/latest/core/engines.html
        """
        filename = kwargs.pop('filename', None)
        obj = cls(**kwargs)
        obj.filename = filename
        obj.connection_uri = connection_uri
        obj._metadata = sql.MetaData(bind=engine)
        obj._initialize_with_mode(self.mode)
        return obj


    def close(self):
        # is this necessary?
        self.engine.dispose()

    @property
    def metadata(self):
        return self._metadata

    @staticmethod
    def filename_from_dialect(filename, dialect):
        # TODO: I think this might be removed in the future; instead, custom
        # setups can use the from_engine constructor

        # take dialects like "sqlite", etc and return proper file connection
        # URI; would be essentially no-op for regular file opening as with
        # .nc
        uri_root = {
            'sqlite': "sqlite:///{filename}",
        }[dialect]
        return uri_root.format(filename=filename)

    def internal_tables_from_db(self):
        """Obtain mappings of table name to number from database.
        """
        tables = self.metadata.tables['tables']
        with self.engine.connect() as conn:
            res = conn.execute(tables.select())
            table_to_number = {r.name: r.idx for r in res}
        number_to_table = {v: k for (k, v) in table_to_number.items()}
        return table_to_number, number_to_table

    def _add_table_to_tables_list(self, table_name, table_schema, cls_):
        if table_name in ['uuid', 'tables']:
            return
        # note that this return the number of tables in 'tables', which does
        # not include 'uuid' or 'tables'
        # There must be a better way to do this, but this seems to get the
        # job done,
        tables = self.metadata.tables['tables']
        with self.engine.connect() as conn:
            res = conn.execute(tables.select())
            n_tables = len(list(res))

        schema_table = self.metadata.tables['schema']

        module = cls_.__module__
        class_name = cls_.__name__

        with self.engine.connect() as conn:
            conn.execute(tables.insert().values(name=table_name,
                                                idx=n_tables,
                                                module=module,
                                                class_name=class_name))
            conn.execute(schema_table.insert().values(
                table=table_name,
                schema=json.dumps(table_schema)
            ))

        self.table_to_number.update({table_name: n_tables})
        self.number_to_table.update({n_tables: table_name})

    def _load_from_table(self, table_name, idx_list):
        # this is not public API (assumes idx_list, which is reserved by not
        # guaranteed)
        table = self.metadata.tables[table_name]
        results = []
        with self.engine.connect() as conn:
            for block in grouper(idx_list, self.MAX_SQL_ITEMS):
                or_stmt = sql.or_(*(table.c.idx == idx for idx in block))
                sel = table.select(or_stmt)
                results.extend(list(conn.execute(sel)))

        return results

    ### FROM HERE IS THE GENERIC PUBLIC API
    def register_type(self, type_str, backend_type):
        self.known_types[type_str] = backend_type

    def register_schema(self, schema, table_to_class,
                        sql_schema_metadata=None):
        """Register (part of) a schema (create necessary tables in DB)

        Raises
        ------
        TypeError
            if the schema provided has names in the existing schema (may be
            trying to modify existing schema)

        """
        for table_name in schema:
            logger.info("Add schema table " + str(table_name))
            columns = make_columns(table_name, schema, sql_schema_metadata,
                                   self.known_types)
            try:
                table = sql.Table(table_name, self.metadata, *columns)
            except sql.exc.InvalidRequestError:
                raise TypeError("Schema registration problem. Your schema "
                                "may already have tables of the same names.")

            if table_name not in universal_schema:
                self._add_table_to_tables_list(table_name,
                                               schema[table_name],
                                               table_to_class[table_name])

        self.metadata.create_all(self.engine)
        self.schema.update(schema)

    def has_table(self, table_name):
        """Returns whether this table is known to the database.

        Parameters
        ----------
        table_name : str
            the name of the table to search for
        """
        return table_name in self.metadata.tables

    def register_storable_function(self, table_name, result_type):
        """
        Parameters
        ----------
        table_name : Str
            the name for this table; typically the UUID of the storable
            function
        result_type : Str
            string name of the result type
        """
        logger.info("Registering storable function: UUID: %s (%s)" %
                    (table_name, result_type))
        col_type = sql_type[backend_registration_type(result_type)]
        columns = [sql.Column('uuid', sql.String, primary_key=True),
                   sql.Column('value', col_type)]
        try:
            table = sql.Table(table_name, self.metadata, *columns)
        except sql.exc.InvalidRequestError:
            raise TypeError("Schema registration problem. Your schema "
                            "may already have tables of the same names.")

        self.metadata.create_all(self.engine)
        # TODO : do we need to do anything else for this?

    def add_storable_function_results(self, table_name, result_dict):
        """
        Parameters
        ----------
        table_name : Str
            name of table for this storable function (typically the UUID)
        result_dict : Mapping[Str, Any]
            mapping from UUID to result
        """
        # anything in the cache doesn't need to be saved
        known_uuids = self.known_uuids[table_name]
        set_uuids = set(result_dict.keys())
        unknown_uuids = set_uuids - known_uuids

        # now we search the for existing
        found_results = self.load_storable_function_results(
            table_name,
            list(unknown_uuids)
        )

        unknown_uuids -= set(found_results.keys())

        # only store the results that haven't been stored
        results = [{'uuid': uuid, 'value': result_dict[uuid]}
                   for uuid in unknown_uuids]
        table = self.metadata.tables[table_name]
        if results:
            with self.engine.connect() as conn:
                conn.execute(table.insert(), results)

        # update the cache
        self.known_uuids[table_name].update(set_uuids)

    def load_storable_function_results(self, table_name, uuids):
        """Load results for given stored function and input UUIDs.

        Parameters
        ----------
        table_name : Str
            name of table for this storable function (typically the UUID)
        uuids : List[Str]
            list of UUIDs to load the results for

        Returns
        -------
        Dict[Str, Any] :
            mapping of UUID to associated value
        """
        table = self.metadata.tables[table_name]
        results = []
        for uuid_block in tools.block(uuids, self.MAX_SQL_ITEMS):
            # uuid_sel = table.select(
                # sql.exists().where(table.c.uuid.in_(uuid_block))
            # )
            uuid_sel = table.select().where(table.c.uuid.in_(uuid_block))
            with self.engine.connect() as conn:
                res = conn.execute(uuid_sel)
                res = res.fetchall()
            results += res

        logger.debug("Found {} UUIDs".format(len(results)))
        result_dict = {uuid: value for uuid, value in results}
        return result_dict

    def load_storable_function_table(self, table_name):
        return {row['uuid']: row['value']
                for row in self.table_iterator(table_name)}

    def add_tag(self, table_name, name, content):
        table = self.metadata.tables[table_name]

        with self.engine.connect() as conn:
            conn.execute(table.insert(), [{'name': name,
                                           'content': content}])


    def add_to_table(self, table_name, objects):
        """Add a list of objects of a given class

        Parameters
        ----------
        table_name : str
            the name of the table
        objects : list of dict
            dict representation of the objects to be added
        """
        # this will insert objects into the table
        table = self.metadata.tables[table_name]
        table_num = self.table_to_number[table_name]

        # this is if we don't use the UUID in the schema... but doing so
        # would be another option (redundant data, but better sanity checks)

        with self.engine.connect() as conn:
            conn.execute(table.insert(), objects)

        res = []
        uuids = [obj['uuid'] for obj in objects]
        for uuid_block in tools.block(uuids, self.MAX_SQL_ITEMS):
            sel_uuids_idx = sql.select([table.c.uuid, table.c.idx]).\
                    where(table.c.uuid.in_(uuid_block))
            with self.engine.connect() as conn:
                res += list(conn.execute(sel_uuids_idx))

        uuid_to_rows = {r[0]: r[1] for r in res}


        uuid_table = self.metadata.tables['uuid']
        uuid_insert_dicts = [{'uuid': k, 'table': table_num, 'row':v}
                             for (k, v) in uuid_to_rows.items()]

        with self.engine.connect() as conn:
            # here we use executemany for performance
            conn.execute(uuid_table.insert(), uuid_insert_dicts)

    def load_n_rows_from_table(self, table_name, first_row, n_rows):
        idx_list = list(range(first_row, first_row + n_rows))
        return self._load_from_table(table_name, idx_list)

    def load_uuids_table(self, uuids, ignore_missing=False):
        """Loads uuids and info on finding data within the table.

        This can also be used to identify which UUIDs already exist in the
        database (with ignore_missing=True).

        Parameters
        ----------
        uuids : list
            list of UUIDs to look up

        Returns
        -------
        list
            table rows with UUID, table ID, and table row index for each
            desired UUID
        """
        uuid_table = self.metadata.tables['uuid']
        logger.debug("Looking for {} UUIDs".format(len(uuids)))
        results = []
        for uuid_block in tools.block(uuids, self.MAX_SQL_ITEMS):
            logger.debug("New block of {} UUIDs".format(len(uuid_block)))
            uuid_sel = uuid_table.select().\
                    where(uuid_table.c.uuid.in_(uuid_block))
            # uuid_or_stmt = sql.or_(*(uuid_table.c.uuid == uuid
                                     # for uuid in uuids))
            # uuid_sel = uuid_table.select(uuid_or_stmt)
            with self.engine.connect() as conn:
                res = list(conn.execute(uuid_sel))
            if not ignore_missing and len(results) != len(uuids):
                # TODO
                # figure out which UUID is missing, raise error on first found
                pass
            results += res

        logger.debug("Found {} UUIDs".format(len(results)))

        return results

    def load_table_data(self, uuid_table_rows):
        # this pulls out a table the information for the relevant UUIDs
        uuid_table_rows = list(uuid_table_rows)  # iterator to list
        by_table_number = tools.group_by_index(uuid_table_rows, 1)
        by_table_name = {self.number_to_table[k]: v
                         for (k, v) in by_table_number.items()}
        loaded_results = []
        for table in by_table_name:
            idxs = [val[2] for val in by_table_name[table]]
            loaded_results += self._load_from_table(table, idxs)

        if self.debug:
            loaded_uuids = [res.uuid for res in loaded_results]
            input_uuids = [row.uuid for row in uuid_table_rows]
            assert set(input_uuids) == set(loaded_uuids)  # sanity check
        return loaded_results

    def database_schema(self):
        """Reload the schema as stored in the database.

        One of the special tables in the database stores the schema
        information. This reloads that into the standard dict form for a
        schema.

        Returns
        -------
        schema : dict
            schema dictionary
        """
        schema_table = self.metadata.tables['schema']
        sel = schema_table.select()
        with self.engine.connect() as conn:
            schema_rows = list(conn.execute(schema_table.select()))
        schema = {r.table: list(map(tuple, json.loads(r.schema)))
                  for r in schema_rows}
        return schema

    def get_representative(self, table_name):
        """Get a row from the table (as a representative of table data)

        This gets a single row from the given table. This can be useful for
        introspection of table structure, or for providing an example object
        (e.g., to use as a template for identifying results of applying a
        function to that kind of object).

        Parameters
        ----------
        table_name : str
            name of the table

        Returns
        -------
        ??? TODO
        """
        table = self.metadata.tables[table_name]
        with self.engine.connect() as conn:
            results = conn.execute(table.select())
            representative = results.fetchone()
            results.close()
        return representative

    @property
    def table_to_class(self):
        tables_table = self.metadata.tables['tables']
        with self.engine.connect() as conn:
            rows = list(conn.execute(tables_table.select()))
        table_to_class = {}
        for row in rows:
            table = row.name
            cls = import_class(row.module, row.class_name)
            table_to_class[table] = cls
        return table_to_class

    def uuid_row_to_table_name(self, uuid_row):
        """Identify the table name from a row from the UUIDs table
        """
        return self.number_to_table[uuid_row.table]

    def table_iterator(self, table_name):
        """Iterate over l rows in the table
        """
        # TODO: this can probably be done in a more low-level way that
        # doesn't require memory caching everything
        table = self.metadata.tables[table_name]
        with self.engine.connect() as conn:
            results = list(conn.execute(table.select()))
        for row in results:
            yield row

    def table_len(self, table_name):
        table = self.metadata.tables[table_name]
        count_query = sql.select([sql.func.count()]).select_from(table)
        with self.engine.connect() as conn:
            results = conn.execute(count_query)
            count_list = [r for r in results]

        if self.debug:
            assert len(count_list) == 1
            subcount = count_list[0]
            assert len(subcount) == 1

        count = count_list[0][0]
        return count

    def table_get_item(self, table_name, item):
        table = self.metadata.tables[table_name]
        # SQL counts from 1; Python counts from 0
        item_sel = table.select().where(table.c.idx == item + 1)
        with self.engine.connect() as conn:
            results = list(conn.execute(item_sel))

        if self.debug:
            assert len(results) == 1

        return results[0]
