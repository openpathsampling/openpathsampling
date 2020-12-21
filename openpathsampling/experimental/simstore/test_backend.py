from .backend import *
import pytest

def test_extract_metadata():
    sql_meta = {
        'uuid': {'uuid': {'primary_key': True}}
    }
    meta = extract_backend_metadata(sql_meta, 'uuid', 'uuid')
    assert meta == {'primary_key': True}

