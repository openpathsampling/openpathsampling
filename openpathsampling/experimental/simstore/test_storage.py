from .storage import *

import pytest

from .serialization_helpers import get_uuid
from .test_utils import all_objects, UnnamedUUID, MockUUIDObject

class TestStorageTable(object):
    def setup(self):
        pass

    def test_iter(self):
        pytest.skip()

    def test_getitem(self):
        pytest.skip()

    def test_len(self):
        pytest.skip()

    def test_save(self):
        pytest.skip()


class TestPseudoTable(object):
    def setup(self):
        self.objs = {None: UnnamedUUID(normal_attr=10)}
        self.objs.update(all_objects)
        self.pseudo_table = PseudoTable(list(self.objs.values()))

    @pytest.mark.parametrize('name', ['int', None])
    def test_get_by_uuid(self, name):
        obj = self.objs[name]
        uuid = get_uuid(obj)
        expected = self.objs[name]
        assert self.pseudo_table.get_by_uuid(uuid) == expected

    @pytest.mark.parametrize('name', ['int', None])
    def test_getitem(self, name):
        obj = self.objs[name]
        uuid = get_uuid(obj)
        idx = self.pseudo_table.index(obj)
        assert self.pseudo_table[idx] == obj
        if name is not None:
            assert self.pseudo_table[name] == obj

    @pytest.mark.parametrize('name', ['int2', None])
    def test_setitem(self, name):
        value = {'int2': 20, None: 21}[name]
        item = {'int2': MockUUIDObject(name="int2", normal_attr=value),
                None: UnnamedUUID(normal_attr=value)}[name]

        len_table = len(self.pseudo_table)

        self.pseudo_table[len_table-1] = item
        assert item in self.pseudo_table
        assert self.pseudo_table[len_table-1] == item
        if name != None:
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
        if name != None:
            assert self.pseudo_table[name] == item

    def test_delitem(self):
        assert len(self.pseudo_table) == len(all_objects) + 1
        assert self.objs['int'] in self.pseudo_table
        int_idx = self.pseudo_table._sequence.index(self.objs['int'])
        del self.pseudo_table[int_idx]
        assert len(self.pseudo_table) == len(all_objects)
        assert self.objs['int'] not in self.pseudo_table
        with pytest.raises(KeyError):
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
        if name != None:
            assert self.pseudo_table[name] == item
