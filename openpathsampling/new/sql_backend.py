import os
import collections
import sqlalchemy as sql
from storage import universal_schema
from tools import group_by, compare_sets
import tools
import ujson as json

from my_types import parse_ndarray_type

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
    'ndarray': sql.LargeBinary,  #TODO: numpy store/load
    #TODO add more
}

universal_sql_meta = {
    'uuid': {'uuid': {'primary_key': True}},
    'tables': {'name': {'primary_key': True}}
}

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
        self.initialize_with_mode(self.mode)

    def initialize_with_mode(self, mode):
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
        obj = cls.__new__(cls)
        self._metadata = sql.MetaData()
        self.schema = {}
        self.table_to_number = {}
        self.number_to_table = {}
        self.engine = engine
        self.mode = mode
        self.initialize_with_mode(self.mode)


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

    @staticmethod
    # TODO: this is not going to be specific to SQL; refactor
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

    def internal_tables_from_db(self):
        """Obtain mappings of table name to number from database.
        """
        tables = self.metadata.tables['tables']
        with self.engine.connect() as conn:
            res = conn.execute(tables.select())
            table_to_number = {r.name: r.idx for r in res}
        number_to_table = {v: k for (k, v) in table_to_number.items()}
        return table_to_number, number_to_table

    def table_list_is_consistent(self):
        """Test whether the DB, schema, and internal table list agree.
        """
        db_tables = set(self.engine.table_names())
        db_tables_tables = self.metadata.tables['tables']
        schema_tables = set(self.schema.keys())
        internal_tables_1 = set(self.table_to_number.keys())
        internal_tables_2 = set(self.number_to_table.values())
        consistent = (db_table == scheme_tables == internal_tables_1
                      == interal_tables_2)
        return consistent

    def table_inconsistencies(self, table_name):
        """Test whether a given table is consistent with its entries in the
        UUID table.
        """
        tables = self.storage.metadata.tables
        table_num = self.table_to_number[table_name]
        uuid_db = tables['uuid']
        uuid_sel = uuid_db.select().where(uuid_db.c.table == table_num)
        with self.engine.connect() as conn:
            table_data = list(conn.execute(tables[table_name].select()))
            uuid_data = list(conn.execute(uuid_sel))
        table_uuids = set([val.uuid for val in table_data])
        uuid_uuids = set([val.uuid for val in uuid_data])
        (table_only, uuid_only) = compare_sets(table_uuids, uuid_uuids)
        return table_only, uuid_only

    def table_is_consistent(self, table_name):
        """
        Are all UUIDs in this table also in the UUID table, and vice versa?

        More specifically, this compares all UUIDs in this table with all
        the objects pointing to this table in the UUID table.

        Note that this is an expensive process, and should be guaranteed by
        the insertion process. But this can act as a useful check on the
        data. The :meth:`.table_inconsistencies` method gives the detailed
        results of the comparison, for debugging purposes.
        """
        table_only, uuid_only = self.table_inconsistencies(table_name)
        if table_only == set([]) and uuid_only == set([]):
            return True
        else:
            return False

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


    def parse_registration_type(self, type_name):
        ops_type = type_name
        ndarray_info = parse_ndarray_type(type_name)
        if parse_ndarray_type(type_name):
            ops_type = 'ndarray'
        return ops_type

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
            columns = []
            if table_name not in ['uuid', 'tables']:
                columns.append(sql.Column('idx', sql.Integer,
                                          primary_key=True))
                columns.append(sql.Column('uuid', sql.String))
            for col, type_name in schema[table_name]:
                # TODO: more general creation of type name
                col_type = sql_type[self.parse_registration_type(type_name)]
                metadata = self._extract_metadata(sql_schema_metadata,
                                                  table_name, col)
                columns.append(sql.Column(col, col_type, **metadata))

            try:
                logger.info("Add schema table" + str(table_name))
                table = sql.Table(table_name, self.metadata, *columns)
            except sql.exc.InvalidRequestError:
                raise TypeError("Schema registration problem. Your schema "
                                "may already have tables of the same names.")

            #TODO: add schema to schema table
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
        # pop_uuids = [{k: v for (k, v) in obj.items() if k != 'uuid'}
                     # for obj in objects]
        insert_statements = [table.insert().values(**obj)
                             for obj in objects]

        with self.engine.connect() as conn:
            # can't use executemany here because we need the resulting
            # primary key values
            uuid_to_rows = {
                obj['uuid']: conn.execute(ins).inserted_primary_key[0]
                for (obj, ins) in zip(objects, insert_statements)
            }

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
        """
        uuid_table = self.metadata.tables['uuid']
        uuid_or_stmt = sql.or_(*(uuid_table.c.uuid == uuid
                                 for uuid in uuids))
        uuid_sel = uuid_table.select(uuid_or_stmt)
        with self.engine.connect() as conn:
            results = list(conn.execute(uuid_sel))
        if not ignore_missing and len(results) != len(uuids):
            # TODO
            # figure out which UUID is missing, raise error on first found
            pass

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
        schema_table = self.metadata.tables['schema']
        sel = schema_table.select()
        with self.engine.connect() as conn:
            schema_rows = list(conn.execute(schema_table.select()))
        schema = {r.table: map(tuple, json.loads(r.schema))
                  for r in schema_rows}
        return schema
