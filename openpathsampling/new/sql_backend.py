import os
import collections
import sqlalchemy as sql
from storage import universal_schema

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
        number_to_table = {v: k for (k, v) in table_to_number.items}
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
        pop_uuids = [{k: v for (k, v) in obj.items() if k != 'uuid'}
                     for obj in objects]
        insert_statements = [table.insert().values(**obj)
                             for obj in pop_uuids]

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
        pass

    def load_uuids(self, uuids, ignore_missing=False):
        """Loads uuids, rows

        This can also be used to identify which UUIDs already exist in the
        database.
        """
        uuid_table = self.metadata.tables['uuid']
        uuid_or_stmt = sql.or_(*(uuid_table.c.uuid == uuid
                                 for uuid in uuids))
        uuid_sel = uuid_table.select(uuid_or_stmt)
        # TODO: something like this
        # pull `results` from the DB
        results = {}
        if not ignore_missing and len(results) != len(uuid):
            # figure out which UUID is missing, raise error on first found
            pass

    def load_table_data(self, uuids):
        # this pulls out a table the information for the relevant UUIDs
        uuid_table_row = self.load_uuids(uuids)
        # group by table
        # load each table
        # return dict of uuid: dict for row

        pass
