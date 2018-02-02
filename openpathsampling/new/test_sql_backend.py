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

    def test_is_consistent(self):
        pytest.skip()

    def test_load_uuids(self):
        pytest.skip()

    def test_load_n_rows_from_table(self):
        pytest.skip()

    def test_add_to_table(self):
        pytest.skip()

    def test_load_table_data(self):
        pytest.skip()

    def test_persistence(self):
        pytest.skip()

