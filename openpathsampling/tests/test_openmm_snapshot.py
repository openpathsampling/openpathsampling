import openpathsampling.engines.openmm as omm_engine
import openpathsampling as paths
from nose.tools import assert_equal
from test_helpers import data_filename

import openmmtools as omt

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
