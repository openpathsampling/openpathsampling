import pytest

import openpathsampling as paths

from openpathsampling.tests.test_helpers import make_1d_traj
from ..simstore.serialization_helpers import get_uuid

from .snapshots import *

SCHEMA = {
    'statics': [('coordinates', 'ndarray.float32({n_atoms},{n_spatial})'),
                ('box_vectors', 'ndarray.float32({n_spatial},{n_spatial})'),
                ('engine', 'uuid')],
    'snapshot': [('statics', 'lazy')]
}

TOY_SCHEMA = {
    'snapshot': [
        ('velocities', 'ndarray.float32({n_atoms},{n_spatial})'),
        ('coordinates', 'ndarray.float32({n_atoms},{n_spatial})'),
        ('engine', 'uuid')
    ]
}

def test_schema_from_entries():
    class Statics(object):
        schema_entries = [('statics', SCHEMA['statics'])]

    assert schema_from_entries([Statics], lazies=['statics']) == SCHEMA

def test_schema_for_snapshot():
    snapshot = make_1d_traj([0.0])[0]
    assert schema_for_snapshot(snapshot) == TOY_SCHEMA

def test_replace_schema_dimensions():
    descriptor = frozenset([('n_atoms', 1000), ('n_spatial', 3),
                            ('class', 'FooSnapshot')])
    expected = {
        'statics': [('coordinates', 'ndarray.float32(1000,3)'),
                    ('box_vectors', 'ndarray.float32(3,3)'),
                    ('engine', 'uuid')],
        'snapshot': [('statics', 'lazy')]
    }
    assert replace_schema_dimensions(SCHEMA, descriptor) == expected

def test_snapshot_registration_from_db():
    pytest.skip()

def test_snapshot_registration_info():
    snapshot = make_1d_traj([0.0])[0]
    schema, class_info_list = snapshot_registration_info(snapshot, 3)
    expected_schema = {'snapshot3': TOY_SCHEMA['snapshot']}
    assert schema == expected_schema
    assert len(class_info_list) == 1
    info = class_info_list[0]
    assert info.table == 'snapshot3'
    assert info.lookup_result == (get_uuid(snapshot.engine),
                                  snapshot.__class__)
    assert info.cls == snapshot.__class__
