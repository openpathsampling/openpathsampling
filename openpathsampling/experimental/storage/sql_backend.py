import os
import collections
import sqlalchemy as sql
from .storage import universal_schema
from .tools import group_by, compare_sets
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
    #TODO add more
}

universal_sql_meta = {
    'uuid': {'uuid': {'primary_key': True}},
    'tables': {'name': {'primary_key': True}}
}

def make_columns(table_name, schema, sql_schema_metadata):
    columns = []
    if table_name not in ['uuid', 'tables']:
        columns.append(sql.Column('idx', sql.Integer,
                                  primary_key=True))
        columns.append(sql.Column('uuid', sql.String))
    for col, type_name in schema[table_name]:
        col_type = sql_type[backend_registration_type(type_name)]
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

class SQLStorageBackend(object):
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
    echo : bool
        whether the engine to echo SQL commands to stdout (useful for
        debugging)
    """
    def __init__(self, filename, mode='r', sql_dialect=None, echo=False):
        self.filename = filename
        sql_dialect = tools.none_to_default(sql_dialect, 'sqlite')
        self.mode = mode
        self.debug = False
        self.max_query_size = 900

        # override later if mode == 'r' or 'a'
        self.schema = {}
        self.table_to_number = {}
        self.number_to_table = {}

        if self.mode == "w" and os.path.exists(filename):
            # delete existing file; write after
            os.remove(filename)

        # we prevent writes by disallowing write method in read mode;
        # for everything else; just connect to the database
        connection_uri = self.filename_from_dialect(filename, sql_dialect)
        self.engine = sql.create_engine(connection_uri,
                                        echo=echo)
        self._metadata = sql.MetaData(bind=self.engine)
        self._initialize_with_mode(self.mode)

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
            self.register_schema(universal_schema, {})
        elif mode == "r" or mode == "a":
            self.metadata.reflect(self.engine)
            self.schema = self.database_schema()
            self.table_to_number, self.number_to_table = \
                    self.internal_tables_from_db()

    @classmethod
    def from_engine(cls, engine, mode='r'):
        """Constructor allowing user to specify the SQLAlchemy Engine.

        The default constructor doesn't allow all the options to the engine
        that SQLAlchemy can provide. Use this if you want more customization
        of the database engine.
        """
        obj = cls.__new__(cls)
        self._metadata = sql.MetaData()
        self.schema = {}
        self.table_to_number = {}
        self.number_to_table = {}
        self.engine = engine
        self.mode = mode
        self._initialize_with_mode(self.mode)


    def close(self):
        # is this necessary?
        self.engine.dispose()

    @property
    def metadata(self):
        return self._metadata

    @staticmethod
    def filename_from_dialect(filename, dialect):
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
        or_stmt = sql.or_(*(table.c.idx == idx for idx in idx_list))
        sel = table.select(or_stmt)
        with self.engine.connect() as conn:
            results = list(conn.execute(sel))
        return results


    ### FROM HERE IS THE GENERIC PUBLIC API
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
            columns = make_columns(table_name, schema, sql_schema_metadata)
            try:
                table = sql.Table(table_name, self.metadata, *columns)
            except sql.exc.InvalidRequestError:
                raise TypeError("Schema registration problem. Your schema "
                                "may already have tables of the same names.")

            if table_name not in ['uuid', 'tables']:
                self._add_table_to_tables_list(table_name,
                                               schema[table_name],
                                               table_to_class[table_name])

        self.metadata.create_all(self.engine)
        self.schema.update(schema)

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
        for uuid_block in tools.block(uuids, self.max_query_size):
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
        for uuid_block in tools.block(uuids, self.max_query_size):
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
        """Iterate over all rows in the table
        """
        table = self.metadata.tables[table_name]
        with self.engine.connect() as conn:
            results = list(conn.execute(table.select()))
        for row in results:
            yield row
