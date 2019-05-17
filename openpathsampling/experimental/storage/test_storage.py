from .storage import *

import pytest

from .serialization_helpers import get_uuid
from .test_utils import all_objects, UnnamedUUID

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
        unnamed = UnnamedUUID(normal_attr=10)
        self.objs = {None: unnamed}
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

    def test_setitem(self):
        pytest.skip()

    def test_delitem(self):
        pytest.skip()

    def test_len(self):
        assert len(self.pseudo_table) == len(all_objects) + 1

    def test_insert(self):
        pytest.skip()
