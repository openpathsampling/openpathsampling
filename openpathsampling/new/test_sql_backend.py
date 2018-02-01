from storage import *
import pytest

class TestSQLStorageBackend(object):
    def setup(self):
        self._delete_tmp_files()
        self.storage = self._default_backend
        self.schema = {
            'uuid': [('uuid', 'uuid'),
                     ('table', 'int'),
                     ('row', 'int')],
            'samples': [('replica', 'int'),
                        ('ensemble', 'uuid'),
                        ('trajectory', 'uuid')]
        }
        self.storage.register_schema(self.schema)

    def teardown(self):
        self._delete_tmp_files()

    @property
    def _default_backend(self):
        # NOTE: test other SQL backends by subclassing and changing this
        return SQLStorageBackend(":memory:")

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
        assert set(table_names) == set(['uuid', 'samples'])
        assert self._col_names_set('uuid') == set(['uuid', 'table', 'row'])
        assert self._col_names_set('samples') == \
                set(['idx', 'replica', 'ensemble', 'trajectory'])

    def test_register_schema(self):
        new_schema = {
            'snapshot0': [('filename', 'str'),
                          ('index', 'int')]
        }
        self.storage.register_schema(new_schema)
        table_names = self.storage.engine.table_names()
        assert set(table_names) == set(['uuid', 'samples', 'snapshot0'])
        assert self._col_names_set('snapshot0') == set(['idx', 'filename',
                                                        'index'])

    def test_register_schema_modify_fails(self):
        with pytest.raises(TypeError):
            self.storage.register_schema(self.schema)

    def test_add_to_table(self):
        pytest.skip()

    def test_load_table_data(self):
        pytest.skip()

    def test_persistence(self):
        pytest.skip()
