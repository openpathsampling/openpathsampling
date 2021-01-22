from __future__ import absolute_import
from builtins import hex
from builtins import object
from nose.tools import raises

from nose.plugins.skip import SkipTest
from .test_helpers import u

import logging

import numpy as np
import numpy.testing as npt
try:
    import openmmtools as omt
except ImportError:
    omt = None

import openpathsampling as paths


import openpathsampling.engines.toy as toy_engine
import openpathsampling.engines.openmm as omm_engine

quiet_loggers = ["initialization", "ensemble", "netcdfplus.objects",
                 "netcdfplus.netcdfplus", "pathmover", "netcdfplus.base"]
for logger in quiet_loggers:
    logging.getLogger("openpathsampling."+logger).setLevel(logging.CRITICAL)


class TestFeatures(object):
    def test_copy_with_replacement_toy(self):
        # test a toy snapshot
        init_coord = np.array([1.0, 2.0])
        init_vel = np.array([3.0, 4.0])
        rep_vel = np.array([3.0, 1.0])
        toy_snap = toy_engine.Snapshot(
            coordinates=init_coord, velocities=init_vel)
        toy_copy = toy_snap.copy_with_replacement(velocities=rep_vel)

        assert(toy_snap.velocities[1] == 4.0)
        assert(toy_copy.velocities[1] == 1.0)

    def test_copy_with_replacement_openmm(self):
        if not paths.integration_tools.HAS_OPENMM or omt is None:
            raise SkipTest
        # test an openmm snapshot
        sys = omt.testsystems.AlanineDipeptideVacuum()
        omm_snap = omm_engine.snapshot_from_testsystem(sys)

        rep_vel = omm_snap.velocities * 2.0
        omm_copy = omm_snap.copy_with_replacement(
            velocities=rep_vel)

        npt.assert_allclose(omm_copy.velocities, rep_vel)

        assert(hex(id(omm_snap.kinetics)) != hex(id(omm_copy.kinetics)))
        assert(hex(id(omm_snap.statics)) == hex(id(omm_copy.statics)))

        omm_copy = omm_snap.copy_with_replacement(
            coordinates=[1.0, 2.0, 3.0] * u.nanometer,
            box_vectors=[1.0, 1.0, 1.0] * u.nanometer)

        assert(hex(id(omm_snap.kinetics)) == hex(id(omm_copy.kinetics)))
        assert(hex(id(omm_snap.statics)) != hex(id(omm_copy.statics)))

    @raises(TypeError)
    def test_parameter_error(self):
        init_coord = np.array([1.0, 2.0])
        init_vel = np.array([3.0, 4.0])
        toy_snap = toy_engine.Snapshot(
            coordinates=init_coord, velocities=init_vel)
        toy_snap.copy_with_replacement(dummy=0)
