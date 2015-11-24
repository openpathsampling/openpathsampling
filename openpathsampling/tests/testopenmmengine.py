'''
@author David W.H. Swenson
'''
import time

from nose.tools import (assert_equal, assert_items_equal)
from nose.plugins.skip import SkipTest

from test_helpers import (true_func, data_filename,
                          assert_equal_array_array,
                          assert_not_equal_array_array)
from openpathsampling.openmm_engine import *
from openpathsampling.snapshot import Snapshot
from openpathsampling.snapshot import Momentum, Configuration

import simtk.openmm as mm
from simtk.openmm import app
from simtk import unit


class testOpenMMEngine(object):
    def setUp(self):
        template = paths.tools.snapshot_from_pdb(data_filename("ala_small_traj.pdb"))
        topology = paths.tools.to_openmm_topology(template)

        # Generated using OpenMM Script Builder
        # http://builder.openmm.org

        forcefield = app.ForceField(
            'amber96.xml',  # solute FF
            'tip3p.xml'     # solvent FF
        )

        # OpenMM System
        system = forcefield.createSystem(
            topology,
            nonbondedMethod=app.PME,
            nonbondedCutoff=1.0*unit.nanometers,
            constraints=app.HBonds,
            rigidWater=True,
            ewaldErrorTolerance=0.0005
        )

        # OpenMM Integrator
        integrator = mm.LangevinIntegrator(
            300*unit.kelvin,
            1.0/unit.picoseconds,
            2.0*unit.femtoseconds
        )
        integrator.setConstraintTolerance(0.00001)

        # Engine options
        options = {
            'nsteps_per_frame': 10,
            'platform': 'fastest',
            'solute_indices' : range(22),
            'n_frames_max' : 5,
            'timestep': 2.0*unit.femtoseconds
        }

        self.engine = paths.OpenMMEngine(
            template,
            system,
            integrator,
            options
        )
        self.engine.initialize()

        context = self.engine.simulation.context
        zero_array = np.zeros((self.engine.n_atoms, 3))
        context.setPositions(self.engine.template.coordinates)
        context.setVelocities(u.Quantity(zero_array, u.nanometers / u.picoseconds))

    def teardown(self):
        pass

    def test_sanity(self):
        pass

    def test_snapshot_get(self):
        snap = self.engine.current_snapshot
        state = self.engine.simulation.context.getState(getVelocities=True,
                                                        getPositions=True)
        pos = state.getPositions(asNumpy=True) / u.nanometers
        vel = state.getVelocities(asNumpy=True) / (u.nanometers / u.picoseconds)
        assert_equal_array_array(snap.coordinates / u.nanometers, pos)
        assert_equal_array_array(snap.velocities / (u.nanometers / u.picoseconds),
                                 vel)

    def test_snapshot_set(self):
        pdb_pos = (self.engine.template.coordinates / u.nanometers)
        testvel = []
        testpos = []
        for i in range(len(pdb_pos)):
            testpos.append(list(np.array(pdb_pos[i]) + 
                                np.array([1.0, 1.0, 1.0]))
                          )
            testvel.append([0.1*i, 0.1*i, 0.1*i])

        self.engine.current_snapshot = Snapshot(
            coordinates=testpos,
            velocities=testvel
        )
        state = self.engine.simulation.context.getState(getPositions=True,
                                                        getVelocities=True)
        sim_coords = state.getPositions(asNumpy=True) / u.nanometers
        sim_vels = state.getVelocities(asNumpy=True) / (u.nanometers/u.picoseconds)

        np.testing.assert_almost_equal(testpos, sim_coords, decimal=5)
        np.testing.assert_almost_equal(testvel, sim_vels, decimal=5)

    def test_generate_next_frame(self):
        snap0 = Snapshot(
            configuration=self.engine.current_snapshot.configuration.copy(),
            momentum=self.engine.current_snapshot.momentum.copy()
        )
        new_snap = self.engine.generate_next_frame()
        old_pos = snap0.coordinates / u.nanometers
        new_pos = new_snap.coordinates / u.nanometers
        old_vel = snap0.velocities / (u.nanometers / u.picoseconds)
        new_vel = new_snap.velocities / (u.nanometers / u.picoseconds)
        assert_equal(old_pos.shape, new_pos.shape)
        assert_equal(old_vel.shape, new_vel.shape)
        assert_not_equal_array_array(old_pos, new_pos)
        assert_not_equal_array_array(old_vel, new_vel)

    def test_generate(self):
        traj = self.engine.generate(self.engine.current_snapshot, [true_func])
        assert_equal(len(traj), self.engine.n_frames_max)

    def test_snapshot_timestep(self):
        assert_equal(self.engine.snapshot_timestep, 20 * u.femtoseconds)

    def test_momentum_setter(self):
        raise SkipTest()
        testvel = []
        for i in range(self.engine.n_atoms):
            testvel.append([0.1*i, 0.1*i, 0.1*i])
        self.engine.momentum = Momentum(velocities=testvel,
                                        kinetic_energy=None)
        np.testing.assert_almost_equal(self.engine.current_snapshot.velocities /
                                       (u.nanometers / u.picoseconds), testvel, decimal=5)

    def test_momentum_getter(self):
        momentum = self.engine.momentum
        state = self.engine.simulation.context.getState(getVelocities=True)
        velocities = state.getVelocities(asNumpy=True)
        assert_equal_array_array(
            momentum.velocities / (u.nanometers / u.picoseconds),
            velocities / (u.nanometers / u.picoseconds)
        )

    def test_configuration_setter(self):
        raise SkipTest()
        pdb_pos = (self.engine.template.coordinates / u.nanometers)
        testpos = []
        for i in range(len(pdb_pos)):
            testpos.append(list(np.array(pdb_pos[i]) + 
                                np.array([1.0, 1.0, 1.0]))
                          )
        self.engine.configuration = Configuration(coordinates=testpos)
        np.testing.assert_almost_equal(self.engine.current_snapshot.coordinates /
                                 u.nanometers, testpos, decimal=5)

    def test_configuration_getter(self):
        config = self.engine.configuration
        state = self.engine.simulation.context.getState(getPositions=True)
        positions = state.getPositions(asNumpy=True)
        assert_equal_array_array(
            config.coordinates / u.nanometers,
            positions / u.nanometers
        )
