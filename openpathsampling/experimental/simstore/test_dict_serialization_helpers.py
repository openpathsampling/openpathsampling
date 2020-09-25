import pytest

from .dict_serialization_helpers import *

class TestTupleKeysSerializers(object):
    class ExampleObj(object):
        def __init__(self, foo, bar):
            self.foo = foo
            self.bar = bar

        def __eq__(self, other):
            return self.foo == other.foo and self.bar == other.bar

        def to_dict(self):
            return {'foo': self.foo, 'bar': self.bar}

        @classmethod
        def from_dict(cls, dct):
            return cls(**dct)

    def setup(self):
        self.foo = {('a', 'b'): [1, 2], ('c', 'd'): [3, 4]}
        self.bar = 3
        self.obj = self.ExampleObj(self.foo, self.bar)
        self.dct = {'foo': {'foo_tuple_keys': [('a', 'b'), ('c', 'd')],
                            'foo_values': [[1, 2], [3, 4]]},
                    'bar': 3}

    def test_tuple_keys_to_dict(self):
        decorated = tuple_keys_to_dict(self.ExampleObj.to_dict, 'foo')
        assert decorated(self.obj) == self.dct

    def test_tuple_keys_from_dict(self):
        decorated = tuple_keys_from_dict(self.ExampleObj.from_dict, 'foo')
        # requires explicit cls because we're not binding to the class!
        assert decorated(self.ExampleObj, self.dct) == self.obj
