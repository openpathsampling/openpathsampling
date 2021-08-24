import pytest

from openpathsampling.integration_tools import _chain_import

@pytest.mark.parametrize('modules', [
    ('foo', 'os.path'),
    ('os.path', 'foo')
])
def test_chain_import(modules):
    import os.path
    mod = _chain_import(*modules)
    assert mod is os.path

def test_chain_import_error():
    with pytest.raises(ImportError, match="bar"):
        _chain_import('foo', 'bar')
