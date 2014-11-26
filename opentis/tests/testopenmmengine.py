'''
@author David W.H. Swenson
'''
import os
from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        assert_almost_equal, raises)
from nose.plugins.skip import Skip, SkipTest
from test_helpers import (true_func, data_filename,
                          assert_equal_array_array,
                          assert_not_equal_array_array)

from opentis.openmm_engine import *
from opentis.snapshot import Snapshot
from opentis.snapshot import Momentum, Configuration

from simtk.unit import femtoseconds, picoseconds, nanometers, kelvin, dalton
from simtk.unit import Quantity
import time

class testOpenMMEngine(object):
    def setUp(self):
        options = {'temperature' : 300.0 * kelvin,
                   'collision_rate' : 1.0 / picoseconds,
                   'timestep' : 2.0 * femtoseconds,
                   'nsteps_per_frame' : 10,
                   'n_frames_max' : 5,
                   'start_time' : time.time(),
                   'fn_initial_pdb' : data_filename("ala_small_traj.pdb"),
                   'platform' : 'fastest',
                   'solute_indices' : range(22), 
                   'forcefield_solute' : 'amber96.xml',
                   'forcefield_solvent' : 'tip3p.xml'
                  }
        self.engine = OpenMMEngine(
            filename=data_filename("openmmengine_test.nc"), 
            topology_file=data_filename("ala_small_traj.pdb"), 
            options=options
        )

        context = self.engine.simulation.context
        zero_array = np.zeros((self.engine.n_atoms, 3))
        context.setPositions(self.engine.template.coordinates)
        context.setVelocities(Quantity(zero_array, nanometers / picoseconds))

    def teardown(self):
        if os.path.isfile(data_filename("openmmengine_test.nc")):
            os.remove(data_filename("openmmengine_test.nc"))

    def test_sanity(self):
        assert_equal(os.path.isfile(data_filename("openmmengine_test.nc")),
                     True)

    def test_equilibrate(self):
        snap0 = Snapshot(
            configuration=self.engine.current_snapshot.configuration.copy(),
            momentum=self.engine.current_snapshot.momentum.copy()
        )
        self.engine.equilibrate(5)
        newsnap = self.engine.current_snapshot
        old_pos = snap0.coordinates
        new_pos = newsnap.coordinates
        for atom_i in range(self.engine.n_atoms):
            # NOTE: we use the engine.storage.topology instead of
            # engine.topology because the latter is a mdtraj, former is
            # openmm; we might change that so that engine.simulation is the
            # only place you find an openmm topol
            res_i = self.engine.storage.topology.atom(atom_i).residue
            if res_i in self.engine.solute_indices:
                assert_items_equal(old_pos[atom_i], new_pos[atom_i])
            else:
                exist_diff = False
                for (old, new) in zip(old_pos[atom_i], new_pos[atom_i]):
                    if old != new:
                        exist_diff = True
                assert_equal(exist_diff, True)

    def test_snapshot_get(self):
        snap = self.engine.current_snapshot
        state = self.engine.simulation.context.getState(getVelocities=True,
                                                        getPositions=True)
        pos = state.getPositions(asNumpy=True) / nanometers
        vel = state.getVelocities(asNumpy=True) / (nanometers / picoseconds)
        assert_equal_array_array(snap.coordinates / nanometers, pos)
        assert_equal_array_array(snap.velocities / (nanometers / picoseconds), 
                                 vel)

    def test_snapshot_set(self):
        pdb_pos = (self.engine.template.coordinates / nanometers)
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
        sim_coords = state.getPositions(asNumpy=True) / nanometers
        sim_vels = state.getVelocities(asNumpy=True) / (nanometers/picoseconds)

        np.testing.assert_almost_equal(testpos, sim_coords, decimal=5)
        np.testing.assert_almost_equal(testvel, sim_vels, decimal=5)

    def test_generate_next_frame(self):
        snap0 = Snapshot(
            configuration=self.engine.current_snapshot.configuration.copy(),
            momentum=self.engine.current_snapshot.momentum.copy()
        )
        new_snap = self.engine.generate_next_frame()
        old_pos = snap0.coordinates / nanometers
        new_pos = new_snap.coordinates / nanometers
        old_vel = snap0.velocities / (nanometers / picoseconds)
        new_vel = new_snap.velocities / (nanometers / picoseconds)
        assert_equal(old_pos.shape, new_pos.shape)
        assert_equal(old_vel.shape, new_vel.shape)
        assert_not_equal_array_array(old_pos, new_pos)
        assert_not_equal_array_array(old_vel, new_vel)

    def test_generate(self):
        self.engine.initialized = True
        traj = self.engine.generate(self.engine.current_snapshot, [true_func])
        assert_equal(len(traj), self.engine.n_frames_max)

    def test_momentum_setter(self):
        testvel = []
        for i in range(self.engine.n_atoms):
            testvel.append([0.1*i, 0.1*i, 0.1*i])
        self.engine.momentum = Momentum(velocities=testvel,
                                        kinetic_energy=None)
        np.testing.assert_almost_equal(self.engine.current_snapshot.velocities /
                                       (nanometers / picoseconds), testvel, decimal=5)

    def test_momentum_getter(self):
        momentum = self.engine.momentum
        state = self.engine.simulation.context.getState(getVelocities=True)
        velocities = state.getVelocities(asNumpy=True)
        assert_equal_array_array(
            momentum.velocities / (nanometers / picoseconds),
            velocities / (nanometers / picoseconds)
        )

    def test_configuration_setter(self):
        pdb_pos = (self.engine.template.coordinates / nanometers)
        testpos = []
        for i in range(len(pdb_pos)):
            testpos.append(list(np.array(pdb_pos[i]) + 
                                np.array([1.0, 1.0, 1.0]))
                          )
        self.engine.configuration = Configuration(coordinates=testpos)
        np.testing.assert_almost_equal(self.engine.current_snapshot.coordinates /
                                 nanometers, testpos, decimal=5)

    def test_configuration_getter(self):
        config = self.engine.configuration
        state = self.engine.simulation.context.getState(getPositions=True)
        positions = state.getPositions(asNumpy=True)
        assert_equal_array_array(
            config.coordinates / nanometers,
            positions / nanometers
        )
