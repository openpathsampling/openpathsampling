from storage import *
import pytest

# NOTE: you can add other SQL backends by subclassing and changing setup
class TestSQLStorageBackend(object):
    def setup(self):
        self._delete_tmp_files()
        self.storage = SQLStorageBackend(":memory:")
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

    @staticmethod
    def _delete_tmp_files():
        if os.path.isfile("test.sql"):
            os.remove("test.sql")

    def test_extract_metadata(self):
        pass

    def test_filename_from_backend(self):
        pass

    def _col_names_set(self, table):
        meta = self.storage.metadata
        return set([col.name for col in meta.tables[table].columns])

    def test_setup(self):
        table_names = self.storage.engine.table_names()
        assert set(table_names) == set(['uuid', 'samples'])
        assert self._col_names_set('uuid') == set(['uuid', 'table', 'row'])
        assert self._col_names_set('samples') == set(['replica', 'ensemble',
                                                      'trajectory'])

    def test_register_schema(self):
        new_schema = {
            'snapshot0': [('filename', 'str'),
                          ('index', 'int')]
        }
        self.storage.register_schema(new_schema)
        table_names = self.storage.engine.table_names()
        assert set(table_names) == set(['uuid', 'samples', 'snapshot0'])
        assert self._col_names_set('snapshot0') == set(['filename',
                                                        'index'])

    def test_register_schema_modify_fails(self):
        # run register_schema on self.schema again
        pass

    def test_add_to_table(self):
        pass

    def test_load_table_data(self):
        pass

    def test_persistence(self):
        pass
