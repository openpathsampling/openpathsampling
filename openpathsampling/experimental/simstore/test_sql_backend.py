from .sql_backend import *
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
        self.table_to_class = {'samples': tuple,
                               'snapshot0': tuple,
                               'snapshot1': tuple}
        self.default_table_names = {'uuid', 'tables', 'schema', 'metadata',
                                    'tags'}

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
        self.database.register_schema(schema, self.table_to_class)
        sample_dict = self._sample_data_dict()
        self.database.add_to_table('samples', sample_dict)
        return sample_dict

    def _add_snapshot_data(self):
        snapshot_schema = {'snapshot0': self.schema['snapshot0']}
        self.database.register_schema(snapshot_schema, self.table_to_class)
        snap_dicts = [{'filename': 'file.trr', 'index': 100,
                       'uuid': 'snapuuid'}]
        self.database.add_to_table('snapshot0', snap_dicts)
        return snap_dicts


    def teardown(self):
        self._delete_tmp_files()

    @property
    def _default_database(self):
        # NOTE: test other SQL backends by subclassing and changing this
        return SQLStorageBackend(":memory:", mode='w')

    @staticmethod
    def _delete_tmp_files():
        tmp_files = ['test.sql', 'test1.sql', 'test2.sql']
        for f in tmp_files:
            if os.path.isfile(f):
                os.remove(f)

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

    @pytest.mark.parametrize('db', ['memory', 'write', 'append'])
    def test_setup(self, db):
        database = {
            'memory': self.database,
            'write': SQLStorageBackend('test1.sql', mode='w'),
            'append': SQLStorageBackend('test2.sql', mode='a'),
        }[db]
        table_names = database.engine.table_names()
        assert set(table_names) == self.default_table_names
        assert self._col_names_set('uuid') == {'uuid', 'table', 'row'}
        assert self._col_names_set('tables') == {'name', 'idx', 'module',
                                                 'class_name'}
        assert self._col_names_set('schema') == {'table', 'schema'}

    def test_register_schema(self):
        new_schema = {
            'snapshot1': [('filename', 'str'), ('index', 'int')]
        }
        self.database.register_schema(new_schema, self.table_to_class)
        table_names = self.database.engine.table_names()
        assert set(table_names) == self.default_table_names | {'snapshot1'}
        assert self._col_names_set('snapshot1') == {'idx', 'uuid',
                                                    'filename', 'index'}

    def test_register_schema_modify_fails(self):
        self.database.register_schema(self.schema, self.table_to_class)
        with pytest.raises(TypeError):
            self.database.register_schema(self.schema, self.table_to_class)

    def test_internal_tables_from_db(self):
        self.database.register_schema(self.schema, self.table_to_class)
        tab2num, num2tab = self.database.internal_tables_from_db()
        tables_db = self.database.metadata.tables['tables']
        with self.database.engine.connect() as conn:
            res = list(conn.execute(tables_db.select()))

        assert len(res) == 2
        num2name = {r.idx: r.name for r in res}
        name2num = {r.name: r.idx for r in res}
        assert tab2num == name2num
        assert num2tab == num2name

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
        self.database.register_schema(schema, self.table_to_class)
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

    def test_load_save_large_numbers(self):
        # mainly a smoke test to ensure that we reload everything correctly
        n_objs = 1001
        snapshot_schema = {'snapshot0': self.schema['snapshot0']}
        self.database.register_schema(snapshot_schema, self.table_to_class)
        assert self.database.table_len('snapshot0') == 0
        snap_dicts = [
            {'filename': 'foo.trr', 'index': idx, 'uuid': 'uuid' + str(idx)}
            for idx in range(n_objs)
        ]
        self.database.add_to_table('snapshot0', snap_dicts)
        assert self.database.table_len('snapshot0') == n_objs
        uuids = [dct['uuid'] for dct in snap_dicts]
        uuid_table_rows = self.database.load_uuids_table(uuids)
        assert len(uuid_table_rows) == n_objs
        reloaded = self.database.load_table_data(uuid_table_rows)
        assert len(reloaded) == n_objs

    def test_load_table_data_missing(self):
        pytest.skip()

    def test_load_table_data_missing_ignored(self):
        pytest.skip()

    def test_persistence(self):
        pytest.skip()

    def test_database_schema(self):
        self.database.register_schema(self.schema, self.table_to_class)
        db_schema = self.database.database_schema()
        assert db_schema == self.schema

    def test_get_representative(self):
        samps = self._add_sample_data()
        snaps = self._add_snapshot_data()
        samp_rep = self.database.get_representative('samples')
        samp_dct = [s for s in samps if s['uuid'] == samp_rep.uuid][0]
        assert samp_rep.replica == samp_dct['replica']
        assert samp_rep.trajectory == samp_dct['trajectory']
        assert samp_rep.ensemble == samp_dct['ensemble']
        snap_rep = self.database.get_representative('snapshot0')
        assert snap_rep.filename == 'file.trr'
        assert snap_rep.index == 100

    def test_table_to_class(self):
        pytest.skip()

    @pytest.mark.parametrize('table', ['samples', 'snapshot0'])
    def test_uuid_row_to_table_name(self, table):
        samps = self._add_sample_data()
        snaps = self._add_snapshot_data()
        input_dict = {'samples': samps, 'snapshot0': snaps}[table]
        uuid_rows = sum([self.database.load_uuids_table([s['uuid']])
                         for s in input_dict], [])
        for row in uuid_rows:
            assert self.database.uuid_row_to_table_name(row) == table

    @pytest.mark.parametrize('table', ['samples', 'snapshot0'])
    def test_table_iterator(self, table):
        samps = self._add_sample_data()
        snaps = self._add_snapshot_data()
        input_dicts = {'samples': samps, 'snapshot0': snaps}[table]
        table_iter = list(self.database.table_iterator(table))

        # test the ordering of the iterator results
        for row in table_iter:
            assert table_iter.index(row) == row.idx - 1

        # test correctness
        for row in table_iter:
            dct = [d for d in input_dicts if d['uuid'] == row.uuid][0]
            for attr in dct:
                assert getattr(row, attr) == dct[attr]

    def test_table_len(self):
        schema = {'samples': [('replica', 'int'),
                             ('ensemble', 'uuid'),
                             ('trajectory', 'uuid')]}
        self.database.register_schema(schema, self.table_to_class)
        assert self.database.table_len('samples') == 0
        sample_dict = self._sample_data_dict()
        self.database.add_to_table('samples', sample_dict)
        assert self.database.table_len('samples') == 3


    def test_table_getitem(self):
        sample_dict = self._add_sample_data()
        for num, samp_dct in enumerate(sample_dict):
            row = self.database.table_get_item('samples', num)
            assert row[0] == num + 1
            # not testing row[1], which should be the UUID
            expected = tuple(samp_dct[k]
                             for k in ['replica', 'ensemble', 'trajectory'])
            assert row[2:] == expected
