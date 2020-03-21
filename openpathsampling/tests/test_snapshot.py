"""
@author David W.H. Swenson
"""
from __future__ import absolute_import
import pytest

from builtins import object
from nose.tools import (assert_equal, assert_not_equal, assert_is, raises,
                        assert_true)
from nose.plugins.skip import Skip, SkipTest
from numpy.testing import assert_allclose
from .test_helpers import CallIdentity, raises_with_message_like, assert_close_unit

import numpy as np
import openpathsampling.engines.features as features

from openpathsampling.engines.snapshot import SnapshotFactory
import openpathsampling as paths

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


class TestFeatures(object):
    def test_coordinates(self):
        Snapshot = SnapshotFactory('TestSnapshot', [features.coordinates], 'A simple testing snapshot')
        compate_attribute(Snapshot, 'coordinates', np.array([0.1, 2.0]), identity)

    def test_velocities(self):
        Snapshot = SnapshotFactory('TestSnapshot', [features.velocities], 'A simple testing snapshot')
        compate_attribute(Snapshot, 'velocities', np.array([0.1, 2.0]), minus)

    def test_box_vectors(self):
        Snapshot = SnapshotFactory('TestSnapshot', [features.box_vectors], 'A simple testing snapshot')
        compate_attribute(Snapshot, 'box_vectors', np.array([0.1, 2.0]), identity)

class TestSnapshotCopy(object):
    def test_copy_none(self):
        if not paths.integration_tools.HAS_OPENMM:
            pytest.skip()
        import openpathsampling.engines.openmm as paths_omm
        # let box_vectors and topology to default to None
        snap = paths_omm.MDSnapshot(coordinates=np.array([[0.0, 0.0, 0.0],
                                                          [0.0, 0.0, 0.0]]),
                                    velocities=np.array([[0.0, 0.0, 0.0],
                                                         [0.0, 0.0, 0.0]]))
        new_snap = snap.copy()
        assert_true(new_snap is not snap)
        assert_true(new_snap.coordinates is not snap.coordinates)
        assert_allclose(new_snap.coordinates, snap.coordinates)
        assert_true(new_snap.box_vectors is snap.box_vectors)
        assert_true(new_snap.box_vectors is None)
        assert_true(new_snap.engine is snap.engine)
