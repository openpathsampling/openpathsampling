import openpathsampling.engines.openmm as omm_engine
import openpathsampling as paths
from nose.tools import assert_equal, assert_almost_equal, assert_not_equal
from test_helpers import data_filename, assert_close_unit

import openmmtools as omt
import simtk.unit as u
import numpy as np

import logging

logging.getLogger('openpathsampling.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.ensemble').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.netcdfplus').setLevel(logging.CRITICAL)

class testOpenMMSnapshot(object):
    def setup(self):
        test_system = omt.testsystems.AlanineDipeptideVacuum()
        self.template = omm_engine.snapshot_from_testsystem(test_system)
        self.engine = omm_engine.Engine(
            topology=self.template.topology,
            system=test_system.system,
            integrator=omt.integrators.VVVRIntegrator()
        )
        self.n_atoms = self.engine.topology.n_atoms
        self.engine.current_snapshot = self.template

    def test_masses_from_file(self):
        masses = self.template.masses
        assert_equal(len(masses), self.n_atoms)

    def test_masses_from_simulation(self):
        sim_snap = self.engine.current_snapshot
        masses = sim_snap.masses
        assert_equal(len(masses), self.n_atoms)

    def test_n_degrees_of_freedom(self):
        assert_equal(self.engine.current_snapshot.n_degrees_of_freedom, 51)

    def test_instantaneous_temperature(self):
        vel_unit = u.nanometers / u.picoseconds
        new_velocities = [[1.0, 0.0, 0.0]] * self.n_atoms * vel_unit
        self.engine.current_snapshot = omm_engine.Snapshot.construct(
            coordinates=self.engine.current_snapshot.coordinates,
            box_vectors=self.engine.current_snapshot.box_vectors,
            velocities=new_velocities,
            engine=self.engine
        )
        test_snap = self.engine.current_snapshot
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
