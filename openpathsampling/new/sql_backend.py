import os
import collections
import sqlalchemy as sql
from storage import universal_schema
from tools import group_by

# dict to convert from OPS string type descriptors to SQL types
sql_type = {
    'uuid': sql.String,
    'list_of_uuid': sql.String,
    'str': sql.String,
    'json': sql.String,
    'int': sql.Integer,
    #TODO add more
}

universal_sql_meta = {
    'uuid': {'uuid': {'primary_key': True}},
    'tables': {'name': {'primary_key': True}}
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

        # override later if mode == 'r' or 'a'
        self._metadata = sql.MetaData()
        self.schema = {}  # TODO: how to load these correctly?
        self.table_to_number = {}
        self.number_to_table = {}

        if backend is None:
            backend = 'sqlite'

        if self.mode == "w" and os.path.exists(filename):
            # delete existing file; write after
            os.remove(filename)

        # we prevent writes by disallowing write method in read mode;
        # for everything else; just connect to the database
        connection_uri = self.filename_from_backend(filename, backend)
        self.engine = sql.create_engine(connection_uri)
        if self.mode == "w":
            self.register_schema(universal_schema)

    @property
    def metadata(self):
        return self._metadata

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

    def is_consistent(self):
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

    def _add_table_to_tables_list(self, table_name):
        if table_name in ['uuid', 'tables']:
            return
        # note that this return the number of tables in 'tables', which does
        # not include 'uuid' or 'tables'
        # There must be a better way to do this, but this seems to get the
        # job done,
        tables = self.metadata.tables['tables']
        with self.engine.connect() as conn:
            res = conn.execute(tables.select())
            n_tables = len([r for r in res])

        with self.engine.connect() as conn:
            conn.execute(tables.insert().values(name=table_name,
                                                idx=n_tables))

        self.table_to_number.update({table_name: n_tables})
        self.number_to_table.update({n_tables: table_name})


    ### FROM HERE IS THE GENERIC PUBLIC API
    def register_schema(self, schema, sql_schema_metadata=None):
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
            columns += [
                sql.Column(
                    col, sql_type[type_name],
                    **self._extract_metadata(sql_schema_metadata,
                                             table_name, col)
                )
                for (col, type_name) in schema[table_name]
            ]
            try:
                table = sql.Table(table_name, self.metadata, *columns)
            except sql.exc.InvalidRequestError:
                raise TypeError("Schema registration problem. Your schema "
                                "may already have tables of the same names.")

            self._add_table_to_tables_list(table_name)

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
        # TODO: I think the sanity checks will be worth it
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

    def _load_from_table(self, table_name, idx_list):
        # this is not public API (assumes idx_list, which is reserved by not
        # guaranteed)
        table = self.metadata.tables[table_name]
        or_stmt = sql.or_(*(table.c.idx == idx for idx in idx_list))
        sel = table.select(or_stmt)
        with self.engine.connect() as conn:
            results = list(conn.execute(sel))
        return results

    def load_n_rows_from_table(self, table_name, first_row, n_rows):
        idx_list = list(range(first_row, first_row + n_rows))
        return self._load_from_table(table_name, idx_list)

    def load_uuids(self, uuids, ignore_missing=False):
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

    def load_table_data(self, uuids):
        # this pulls out a table the information for the relevant UUIDs
        uuid_table_row = self.load_uuids(uuids)
        by_table_number = group_by(uuid_table_row, 1)
        by_table_name = {self.number_to_table[k]: v
                         for (k, v) in by_table_number.items()}
        loaded_results = []
        for table in by_table_name:
            idxs = [val[2] for val in by_table_name[table]]
            loaded_results += self._load_from_table(table, idxs)

        loaded_uuids = [res.uuid for res in loaded_results]
        assert set(uuids) == set(loaded_uuids)  # sanity check
        return loaded_results
