from .storage import *

import pytest

from .serialization_helpers import get_uuid
from .test_utils import all_objects, UnnamedUUID, MockUUIDObject, MockStorage


class TestStorageTable(object):
    def setup(self):
        self.storage = MockStorage()
        # Override default tables/uuids
        self.storage.backend.table_names = ["all"]
        self.obj_list = list(all_objects.values())
        self.storage.backend.tables = [self.obj_list]
        self.storage.backend.uuid_table = {
            e.uuid:
            self.storage.backend.row_types["uuids"](uuid=e.uuid, table=0,
                                                    idx=i)
            for i, e in enumerate(self.obj_list)}
        self.table = StorageTable(self.storage, "all")

    def test_iter(self):
        pytest.skip()

    def test_getitem(self):
        pytest.skip()

    def test_len(self):
        pytest.skip()

    def test_save(self):
        pytest.skip()

    @pytest.mark.parametrize("s", [slice(None),  # [:]
                                   slice(5),  # [:5]
                                   slice(-5),  # [:-5]
                                   slice(1, 2),  # [1:2]
                                   slice(100),  # longer slice than possible
                                   slice(6, 4, -1),  # [6:4:-1]
                                   slice(5, 5),  # empty
                                   ])
    def test_getslice(self, s):
        truths = self.obj_list[s]
        tests = self.table[s]
        for i, j in zip(truths, tests):
            assert i == j

    def test_bogus_access(self):
        with pytest.raises(TypeError, match="type tuple"):
            self.table[(1, 2, 3)]


class TestPseudoTable(object):
    def setup(self):
        self.objs = {None: UnnamedUUID(normal_attr=10)}
        self.objs.update(all_objects)
        self.obj_list = list(self.objs.values())
        self.pseudo_table = PseudoTable(self.obj_list)

    @pytest.mark.parametrize('name', ['int', None])
    def test_get_by_uuid(self, name):
        obj = self.objs[name]
        uuid = get_uuid(obj)
        expected = self.objs[name]
        assert self.pseudo_table.get_by_uuid(uuid) == expected

    @pytest.mark.parametrize('name', ['int', None])
    def test_getitem(self, name):
        obj = self.objs[name]
        idx = self.pseudo_table.index(obj)
        assert self.pseudo_table[idx] == obj
        if name is not None:
            assert self.pseudo_table[name] == obj

    @pytest.mark.parametrize("s", [slice(None),  # [:]
                                   slice(5),  # [:5]
                                   slice(-5),  # [:-5]
                                   slice(1, 2),  # [1:2]
                                   slice(100),  # longer slice than possible
                                   slice(6, 4, -1),  # [6:4:-1]
                                   slice(5, 5),  # empty
                                   ])
    def test_getslice(self, s):
        truths = self.obj_list[s]
        tests = self.pseudo_table[s]
        for i, j in zip(truths, tests):
            assert i == j

    def test_bogus_access(self):
        with pytest.raises(TypeError, match="tuple"):
            self.pseudo_table[(1, 2, 3)]

    @pytest.mark.parametrize('name', ['int2', None])
    def test_setitem(self, name):
        value = {'int2': 20, None: 21}[name]
        item = {'int2': MockUUIDObject(name="int2", normal_attr=value),
                None: UnnamedUUID(normal_attr=value)}[name]

        len_table = len(self.pseudo_table)

        self.pseudo_table[len_table-1] = item
        assert item in self.pseudo_table
        assert self.pseudo_table[len_table-1] == item
        if name is not None:
            assert self.pseudo_table[name] == item

    @pytest.mark.parametrize('name', ['int2', None])
    def test_append(self, name):
        value = {'int2': 20, None: 21}[name]
        item = {'int2': MockUUIDObject(name="int2", normal_attr=value),
                None: UnnamedUUID(normal_attr=value)}[name]

        len_table = len(self.pseudo_table)

        self.pseudo_table.append(item)
        assert item in self.pseudo_table
        assert self.pseudo_table[len_table] == item
        if name is not None:
            assert self.pseudo_table[name] == item

    def test_delitem(self):
        assert len(self.pseudo_table) == len(all_objects) + 1
        assert self.objs['int'] in self.pseudo_table
        int_idx = self.pseudo_table._sequence.index(self.objs['int'])
        del self.pseudo_table[int_idx]
        assert len(self.pseudo_table) == len(all_objects)
        assert self.objs['int'] not in self.pseudo_table
        with pytest.raises(KeyError, match='int'):
            self.pseudo_table['int']

    def test_len(self):
        assert len(self.pseudo_table) == len(all_objects) + 1

    @pytest.mark.parametrize('name', ['int2', None])
    def test_insert(self, name):
        value = {'int2': 20, None: 21}[name]
        item = {'int2': MockUUIDObject(name="int2", normal_attr=value),
                None: UnnamedUUID(normal_attr=value)}[name]

        len_table = len(self.pseudo_table)

        self.pseudo_table.insert(len_table, item)
        assert item in self.pseudo_table
        assert self.pseudo_table[len_table] == item
        if name is not None:
            assert self.pseudo_table[name] == item
