from sql_backend import *
import pytest

class TestSQLStorageBackend(object):
    def setup(self):
        self._delete_tmp_files()
        self.database = self._default_database
        self.database.debug = True
        self.schema = {
            'samples': [('replica', 'int'),
                        ('ensemble', 'uuid'),
                        ('trajectory', 'uuid')],
            'snapshot0': [('filename', 'str'), ('index', 'int')]
        }
        self.default_table_names = {'uuid', 'tables', 'schema', 'metadata'}

    def _sample_data_dict(self):
        sample_list = [(0, 'ens1', 'traj1'),
                       (1, 'ens2', 'traj2'),
                       (0, 'ens1', 'traj2')]
        sample_dict = [
            {'replica': s[0], 'ensemble': s[1], 'trajectory': s[2]}
            for s in sample_list
        ]
        for s in sample_dict:
            s.update({'uuid': str(hex(hash(str(sample_dict))))})
        return sample_dict

    def _add_sample_data(self):
        # use this as a fixture when you need sample data
        schema = {'samples': [('replica', 'int'),
                             ('ensemble', 'uuid'),
                             ('trajectory', 'uuid')]}
        self.database.register_schema(schema)
        sample_dict = self._sample_data_dict()
        self.database.add_to_table('samples', sample_dict)
        return sample_dict

    def teardown(self):
        self._delete_tmp_files()

    @property
    def _default_database(self):
        # NOTE: test other SQL backends by subclassing and changing this
        return SQLStorageBackend(":memory:", mode='w')

    @staticmethod
    def _delete_tmp_files():
        if os.path.isfile("test.sql"):
            os.remove("test.sql")

    def test_extract_metadata(self):
        sql_meta = {
            'uuid': {'uuid': {'primary_key': True}}
        }
        meta = self.database._extract_metadata(sql_meta, 'uuid', 'uuid')
        assert meta == {'primary_key': True}

    @pytest.mark.parametrize('test_input,expected', [
        (("file.sql", "sqlite"), "sqlite:///file.sql"),
        ((":memory:", "sqlite"), "sqlite:///:memory:")
    ])
    def test_filename_from_dialect(self, test_input, expected):
        filename_from_dialect = self.database.filename_from_dialect
        assert filename_from_dialect(*test_input) == expected

    def test_filename_from_dialect_unknown(self):
        with pytest.raises(KeyError):
            self.database.filename_from_dialect("file.sql", "foo")

    def _col_names_set(self, table):
        meta = self.database.metadata
        return set([col.name for col in meta.tables[table].columns])

    def test_setup(self):
        table_names = self.database.engine.table_names()
        assert set(table_names) == self.default_table_names
        assert self._col_names_set('uuid') == {'uuid', 'table', 'row'}
        assert self._col_names_set('tables') == {'name', 'idx'}
        assert self._col_names_set('schema') == {'table', 'schema'}

    def test_register_schema(self):
        new_schema = {
            'snapshot1': [('filename', 'str'), ('index', 'int')]
        }
        self.database.register_schema(new_schema)
        table_names = self.database.engine.table_names()
        assert set(table_names) == self.default_table_names | {'snapshot1'}
        assert self._col_names_set('snapshot1') == {'idx', 'uuid',
                                                    'filename', 'index'}

    def test_register_schema_modify_fails(self):
        self.database.register_schema(self.schema)
        with pytest.raises(TypeError):
            self.database.register_schema(self.schema)

    def test_internal_tables_from_db(self):
        self.database.register_schema(self.schema)
        tab2num, num2tab = self.database.internal_tables_from_db()
        tables_db = self.database.metadata.tables['tables']
        with self.database.engine.connect() as conn:
            res = list(conn.execute(tables_db.select()))

        assert len(res) == 2
        num2name = {r.idx: r.name for r in res}
        name2num = {r.name: r.idx for r in res}
        assert tab2num == name2num
        assert num2tab == num2name

    def test_table_list_is_consistent(self):
        pytest.skip()

    def test_table_inconsistencies(self):
        pytest.skip()

    def test_table_is_consistent(self):
        pytest.skip()

    def test_load_uuids_table(self):
        sample_dict = self._add_sample_data()
        uuids = [s['uuid'] for s in sample_dict]
        loaded_uuids = self.database.load_uuids_table(uuids)
        assert len(loaded_uuids) == 3
        for uuid in loaded_uuids:
            assert uuid.uuid in uuids
            assert uuid.table == 0

    def test_load_n_rows_from_table(self):
        pytest.skip()

    def test_add_to_table(self):
        schema = {'samples': [('replica', 'int'),
                             ('ensemble', 'uuid'),
                             ('trajectory', 'uuid')]}
        self.database.register_schema(schema)
        sample_dict = self._sample_data_dict()
        self.database.add_to_table('samples', sample_dict)

        # load created data
        tables = self.database.metadata.tables
        with self.database.engine.connect() as conn:
            samples = list(conn.execute(tables['samples'].select()))
            uuids = list(conn.execute(tables['uuid'].select()))

        # tests
        assert len(samples) == 3
        assert len(uuids) == 3
        for uuid in uuids:
            assert uuid.table == 0
        # check that row numbers match to right UUID
        uuid_dict = {u.uuid: u.row for u in uuids}
        samples_uuid_dict = {s.idx: s.uuid for s in samples}
        for uuid_val in uuid_dict:
            assert samples_uuid_dict[uuid_dict[uuid_val]] == uuid_val

        # check that we got back the objects we created
        returned_sample_dict = [
            {'uuid': s.uuid, 'replica': s.replica, 'ensemble': s.ensemble,
             'trajectory': s.trajectory}
            for s in samples
        ]
        # effectively test sets; but dict isn't hashable; length already
        # tested above
        for dct in returned_sample_dict:
            assert dct in sample_dict

    def test_load_table_data(self):
        sample_dict = self._add_sample_data()
        uuids = [s['uuid'] for s in sample_dict]
        uuid_table_rows = self.database.load_uuids_table(uuids)
        loaded_table = self.database.load_table_data(uuid_table_rows)
        assert len(loaded_table) == 3
        pytest.skip()

    def test_load_table_data_missing(self):
        pytest.skip()

    def test_load_table_data_missing_ignored(self):
        pytest.skip()

    def test_load_table_data_lazy(self):
        pytest.skip()

    def test_persistence(self):
        pytest.skip()

    def test_database_schema(self):
        self.database.register_schema(self.schema)
        db_schema = self.database.database_schema()
        assert db_schema == self.schema

