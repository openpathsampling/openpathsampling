from sql_backend import *
import pytest

class TestSQLStorageBackend(object):
    def setup(self):
        self._delete_tmp_files()
        self.storage = self._default_backend
        self.schema = {
            'samples': [('replica', 'int'),
                        ('ensemble', 'uuid'),
                        ('trajectory', 'uuid')],
            'snapshot0': [('filename', 'str'), ('index', 'int')]
        }

    def teardown(self):
        self._delete_tmp_files()

    @property
    def _default_backend(self):
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
        meta = self.storage._extract_metadata(sql_meta, 'uuid', 'uuid')
        assert meta == {'primary_key': True}

    @pytest.mark.parametrize('test_input,expected', [
        (("file.sql", "sqlite"), "sqlite:///file.sql"),
        ((":memory:", "sqlite"), "sqlite:///:memory:")
    ])
    def test_filename_from_backend(self, test_input, expected):
        filename_from_backend = self.storage.filename_from_backend
        assert filename_from_backend(*test_input) == expected

    def test_filename_from_backend_unknown(self):
        with pytest.raises(KeyError):
            self.storage.filename_from_backend("file.sql", "foo")

    def _col_names_set(self, table):
        meta = self.storage.metadata
        return set([col.name for col in meta.tables[table].columns])

    def test_setup(self):
        table_names = self.storage.engine.table_names()
        assert set(table_names) == set(['uuid', 'tables'])
        assert self._col_names_set('uuid') == set(['uuid', 'table', 'row'])
        assert self._col_names_set('tables') == set(['name', 'idx'])

    def test_register_schema(self):
        new_schema = {
            'snapshot1': [('filename', 'str'), ('index', 'int')]
        }
        self.storage.register_schema(new_schema)
        table_names = self.storage.engine.table_names()
        assert set(table_names) == set(['uuid', 'tables', 'snapshot1'])
        assert self._col_names_set('snapshot1') == set(['idx', 'uuid',
                                                        'filename',
                                                        'index'])

    def test_register_schema_modify_fails(self):
        self.storage.register_schema(self.schema)
        with pytest.raises(TypeError):
            self.storage.register_schema(self.schema)

    def test_internal_tables_from_db(self):
        self.storage.register_schema(self.schema)
        tab2num, num2tab = self.storage.internal_tables_from_db()
        tables_db = self.storage.metadata.tables['tables']
        with self.storage.engine.connect() as conn:
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
        pytest.skip()

    def test_load_n_rows_from_table(self):
        pytest.skip()

    def test_add_to_table(self):
        schema = {'samples': [('replica', 'int'),
                             ('ensemble', 'uuid'),
                             ('trajectory', 'uuid')]}
        self.storage.register_schema(schema)
        sample_list = [(0, 'ens1', 'traj1'),
                       (1, 'ens2', 'traj2'),
                       (0, 'ens1', 'traj2')]
        sample_dict = [
            {'replica': s[0], 'ensemble': s[1], 'trajectory': s[2]}
            for s in sample_list
        ]
        for s in sample_dict:
            s.update({'uuid': str(hex(hash(str(sample_dict))))})
        self.storage.add_to_table('samples', sample_dict)

        # load created data
        tables = self.storage.metadata.tables
        with self.storage.engine.connect() as conn:
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

        pytest.skip()

    def test_load_table_data(self):
        pytest.skip()

    def test_persistence(self):
        pytest.skip()

