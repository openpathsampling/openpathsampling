import pytest

from .write_code import *
from .write_code import (_make_arguments, _call_super, _set_attribute,
                         _def_function)

from openpathsampling.experimental.attach_features import Parameter

@pytest.fixture
def parameters():
    foo = Parameter('foo')
    bar = Parameter('bar', default=None)
    baz = Parameter('baz', default=10)
    qux = Parameter('qux', default='quux')
    return [foo, bar, baz, qux]

def test_make_arguments(parameters):
    expected = "foo, bar=None, baz=10, qux='quux'"
    assert _make_arguments(parameters) == expected

@pytest.mark.parametrize('assign_to', [None, 'foo'])
@pytest.mark.parametrize('params', [True, False])
def test_call_super(assign_to, params, parameters):
    params = parameters if params else []
    param_str = "foo=foo, bar=bar, baz=baz, qux=qux" if params else ""
    assign_str = f"{assign_to} = " if assign_to is not None else ""
    expected = f"{assign_str}super().method_name({param_str})"
    assert _call_super('method_name', params, assign_to) == expected

def test_set_attribute():
    assert _set_attribute('foo') == "self.foo = foo"

@pytest.mark.parametrize('params', [True, False])
def test_make_init_code(params):
    pytest.skip()

def test_make_copy_with_replacement_code():
    pytest.skip()

def test_make_init():
    pytest.skip()

def test_make_copy_with_replacement():
    pytest.skip()
