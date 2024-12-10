import pytest

from .core import *

import types


class TestParameter:
    def setup(self):
        self.param_all_defaults = Parameter('foo')
        self.param_default_none = Parameter('foo', default=None)
        self.param_default_int = Parameter('foo', default=10)
        self.param_default_float = Parameter('foo', default=1.5)
        self.param_default_string = Parameter('foo', default='bar')

    @pytest.mark.parametrize('replace', [
        'name', 'default', 'lazy', 'copy', 'operations', 'docstring'
    ])
    def test_copy_with_replacement(self, replace):
        old = self.param_all_defaults
        replacements = {
            'name': 'bar',
            'default': 'bar',
            'lazy': True,
            'copy': lambda x: 'bar',
            'operations': {'time_reverse': lambda x: 'bar'},
            'docstring': None
        }
        replacement = replacements[replace]
        kwargs = {replace: replacement}
        new = old.copy_with_replacement(**kwargs)
        assert new is not old
        for attr in replacements:
            new_val =  getattr(new, attr)
            if attr == replace:
                assert new_val == replacement
            else:
                assert new_val == getattr(old, attr)

    @pytest.mark.parametrize('default_type', [
        'int', 'float', 'string', 'None', 'no-default'
    ])
    def test_default_as_code(self, default_type):
        default, expected = {
            'int': (10, '10'),
            'float': (1.0, '1.0'),
            'string': ('bar', "'bar'"),
            'None': (None, 'None'),
            'no-default': (Parameter.empty, None),
        }[default_type]
        param = Parameter('foo', default=default)
        assert param.default_as_code == expected

    def test_default_as_code_custom(self):
        import os
        param = Parameter('foo', default=os.getcwd)
        param.default_as_code = 'os.getcwd'
        assert param.default_as_code == 'os.getcwd'

    def test_default_as_code_custom_error(self):
        param = Parameter('foo', default='bar')
        with pytest.raises(TypeError, match="must be a string"):
            param.default_as_code = 10


class TestFunctionFeature:
    @function_feature
    def foo(self, inputs):
        return inputs

    def test_init(self):
        assert isinstance(self.foo, function_feature)
        assert isinstance(self.foo.func, types.FunctionType)
        assert self.foo.func('placeholder', 'bar') == 'bar'

    def test_call(self):
        assert self.foo('placeholder', 'bar') == 'bar'


feature_module_1 = '''
# normally this would require importing Parameter and function_feature, but
# we'll exec in a namespace that already includes those

foo = Parameter('foo',
                docstring=("""
                           foo: int
                               description of foo
                           """))

@function_feature
def snapshot_with_value_added_to_foo(self, value):
    return self.copy_with_replacement(foo=self.foo + value)

@property
def twice_foo(self):
    return 2 * self.foo

# we'll also check that other functions don't get attached
def unattached_function(parameters):
    pass
'''

feature_module_2 = '''
bar = Parameter('bar',
                docstring=("""
                           bar: int
                              description of bar
                           """))
'''

old_style_module = '''
"""
foo : int
    old style docstring
"""

variables = ['foo', 'bar']
numpy = ['foo']  # this affects how copying is performed
minus = ['foo']  # on time reversal, reversed_foo = -foo
flip = ['bar']  # on time reversal, reversed_bar = ~bar (toggle True/False)
lazy = ['foo']  # foo is loaded from storage using a proxy
'''

mixed_style_module = old_style_module + '''
# below is the same behavior as old-style (different docstring)
import numpy as np
import operators

foo = Parameter('foo',
                copy=np.copy,
                operations={'time_reverse': operators.neg},
                lazy=True,
                docstring="""
                          foo : int
                              new style docstring
                          """)

bar = Parameter('bar',
                operations={'time_reverse': operators.invert})
'''

class TestFeatureCollection:
    def test_make_parameters_from_old(self):
        pass

    def test_make_functions_from_old(self):
        pass

    def test_add(self):
        pass

    def test_all_names(self):
        pass

    def test_error_if_overlap_pass(self):
        pass

    def test_error_if_overlap_fail(self):
        pass

    def test_from_module(self):
        pass


def test_attach_features_decorator():
    pass

