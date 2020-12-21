import pytest

from .tools import *

def test_none_to_default():
    assert none_to_default(option=None, default="foo") == "foo"
    assert none_to_default(option="bar", default="foo") == "bar"

def test_compare_sets():
    pytest.skip()

class TestGroupBy(object):
    def setup(self):
        pass

class TestListify(object):
    # TODO: what should be the correct behavior with a dict?
    @pytest.mark.parametrize("obj", [3, "foo"])
    def test_should_wrap(self, obj):
        assert listify(obj) == [obj]

    @pytest.mark.parametrize("obj", [(1, 2), ['a'], {5}])
    def test_should_not_wrap(self, obj):
        assert listify(obj) is obj


class TestFlatten(object):
    def setup(self):
        self.result = ['a', 'b', 'c', 'd', 'e', 'f']
        self.list = ['a', 'b', [['c'], ['d'], 'e'], [['f']]]
        self.dict = {11: 'a', 12: 'b',
                     13: {21: {31: 'c'}, 22: {32: 'd'}, 23: 'e'},
                     14: {24: 'f'}}
        self.mixed_dict = {11: 'a', 12: 'b',
                           13: [{31: 'c'}, {32: 'd'}, 'e'],
                           14: ['f']}
        self.mixed_list = ['a', 'b', {21: ['c'], 22: ['d'], 23: 'e'},
                           {24: 'f'}]

    def test_flatten_dict(self):
        assert flatten_dict(self.dict) == self.result
        expected_mixed = ['a', 'b', [{31: 'c'}, {32: 'd'}, 'e'], ['f']]
        assert flatten_dict(self.mixed_dict) == expected_mixed

    def test_flatten_iterable(self):
        assert flatten_iterable(self.list) == self.result
        assert flatten_iterable(self.mixed_list) == self.mixed_list

    def test_flatten_all(self):
        assert flatten_all(self.list) == self.result
        assert flatten_all(self.dict) == self.result
        assert flatten_all(self.mixed_dict) == self.result
        assert flatten_all(self.mixed_list) == self.result

class TestNestedUpdate(object):
    pass
