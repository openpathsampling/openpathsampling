from .storage import *

import pytest

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
        self.unnamed = UnnamedUUID(normal_attr=10)
        objs = [self.unnamed] + list(all_objects.values())
        self.pseudo_table = PseudoTable(objs)

    def test_get_by_uuid(self):
        pytest.skip()

    def test_getitem(self):
        pytest.skip()

    def test_setitem(self):
        pytest.skip()

    def test_delitem(self):
        pytest.skip()

    def test_len(self):
        pytest.skip()

    def test_insert(self):
        pytest.skip()
