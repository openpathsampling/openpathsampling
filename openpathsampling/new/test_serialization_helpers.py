from serialization_helpers import *
import numpy as np
import pytest

class MockUUIDObject(object):
    def __init__(self, uuid):
        self.__uuid__ = uuid
        self.none_attr = None
        self.dict_attr = {}
        self.list_attr = []
        self.ndarray_attr = np.array([1.0, 1.0])
        self.int_attr = 5
        self.str_attr = "foo"
        self.obj_attr = None
    pass

def test_get_uuid():
    pytest.skip()
