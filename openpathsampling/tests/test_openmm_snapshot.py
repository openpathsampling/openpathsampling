from __future__ import division
from __future__ import absolute_import
from builtins import object
import openpathsampling.engines.openmm as omm_engine
import openpathsampling as paths
from nose.tools import (assert_equal, assert_almost_equal, assert_not_equal,
                        assert_is_not, assert_is)
from nose.plugins.skip import SkipTest
from .test_helpers import data_filename, assert_close_unit, u

import pytest

try:
    import openmmtools as omt
except ImportError:
    omt = None

from openpathsampling.integration_tools import openmm
import numpy as np

import logging

logging.getLogger('openpathsampling.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.ensemble').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.netcdfplus').setLevel(logging.CRITICAL)

class TestOpenMMSnapshot(object):
    def setup(self):
        if not openmm:
            raise SkipTest("OpenMM not installed")
        if not omt:
            raise SkipTest("OpenMMTools not installed; required for OpenMM "
                           "tests.")
        self.test_system = omt.testsystems.AlanineDipeptideVacuum()
        self.template = omm_engine.snapshot_from_testsystem(self.test_system)
        self.engine = omm_engine.Engine(
            topology=self.template.topology,
            system=self.test_system.system,
            integrator=omt.integrators.VVVRIntegrator()
        )
        self.n_atoms = self.engine.topology.n_atoms
        self.test_snap = omm_engine.Snapshot.construct(
            coordinates=self.template.coordinates,
            box_vectors=self.template.box_vectors,
            velocities=self.template.velocities,
            engine=self.engine
        )

    def test_masses_from_file(self):
        masses = self.template.masses
        assert_equal(len(masses), self.n_atoms)

    def test_masses_from_simulation(self):
        assert len(self.test_snap.masses) == self.n_atoms

    def test_n_degrees_of_freedom(self):
        assert self.test_snap.n_degrees_of_freedom == 51

    def test_instantaneous_temperature(self):
        vel_unit = u.nanometers / u.picoseconds
        new_velocities = [[1.0, 0.0, 0.0]] * self.n_atoms * vel_unit
        test_snap = omm_engine.Snapshot.construct(
            coordinates=self.template.coordinates,
            box_vectors=self.template.box_vectors,
            velocities=new_velocities,
            engine=self.engine
        )

        expected_ke = sum(
            [m * vel_unit**2 for m in test_snap.masses],
            0.0*u.joule
        )
        n_dofs = 51.0  # see test above
        assert_close_unit(test_snap.instantaneous_temperature,
                          expected_ke / u.BOLTZMANN_CONSTANT_kB / n_dofs)

    def test_instantaneous_temperature_changes(self):
        trajectory = self.engine.generate(self.template,
                                          [lambda t, foo: len(t) < 4])
        temp_1 = trajectory[1].instantaneous_temperature
        temp_2 = trajectory[2].instantaneous_temperature
        assert_not_equal(temp_1, temp_2)

    def test_mdtraj_trajectory(self):
        snap_1 = omm_engine.snapshot_from_testsystem(self.test_system,
                                                     periodic=False)
        assert_is(snap_1.box_vectors, None)
        traj_1 = snap_1.md
        assert_equal(len(traj_1), 1)
        assert_is_not(traj_1.xyz, None)
        assert_is(traj_1.unitcell_vectors, None)

        snap_2 = self.test_snap
        traj_2 = snap_2.md
        assert_equal(len(traj_2), 1)
        assert_is_not(traj_2.xyz, None)
        assert_is_not(traj_2.unitcell_vectors, None)
