"""
@author David W.H. Swenson
"""

from nose.tools import assert_equal, assert_not_equal, assert_is, raises
from nose.plugins.skip import Skip, SkipTest
from test_helpers import CallIdentity, raises_with_message_like, assert_close_unit

import numpy as np
import openpathsampling.engines.features as features

from openpathsampling.engines.snapshot import SnapshotFactory

def compate_attribute(snapshot_class, attr_name, attr_value, attr_reversal_fnc):
    a = snapshot_class(**{attr_name: attr_value})
    b = a.copy()
    c = a.create_reversed()
    d = a.reversed

    assert_close_unit(getattr(a, attr_name), attr_value)
    assert(getattr(a, attr_name) is not getattr(b, attr_name))
    assert_close_unit(getattr(a, attr_name), getattr(b, attr_name))
    assert_close_unit(getattr(a, attr_name), attr_reversal_fnc(getattr(c, attr_name)))
    assert_close_unit(getattr(a, attr_name), attr_reversal_fnc(getattr(d, attr_name)))

identity = lambda x: x
minus = lambda x: -x
flip = lambda x: not x


class testFeatures(object):
    def test_coordinates(self):
        Snapshot = SnapshotFactory('TestSnapshot', [features.coordinates], 'A simple testing snapshot')
        compate_attribute(Snapshot, 'coordinates', np.array([0.1, 2.0]), identity)

    def test_velocities(self):
        Snapshot = SnapshotFactory('TestSnapshot', [features.velocities], 'A simple testing snapshot')
        compate_attribute(Snapshot, 'velocities', np.array([0.1, 2.0]), minus)

    def test_box_vectors(self):
        Snapshot = SnapshotFactory('TestSnapshot', [features.box_vectors], 'A simple testing snapshot')
        compate_attribute(Snapshot, 'box_vectors', np.array([0.1, 2.0]), identity)
