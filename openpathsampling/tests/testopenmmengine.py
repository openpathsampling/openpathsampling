"""
@author David W.H. Swenson
"""

import numpy as np
import simtk.openmm as mm
from nose.tools import (assert_equal)
from simtk import unit as u
from simtk.openmm import app

import openpathsampling.engines.openmm as peng

from test_helpers import (true_func, data_filename,
                          assert_equal_array_array,
                          assert_not_equal_array_array)


def setUp():
    global topology, template, system
    template = peng.snapshot_from_pdb(data_filename("ala_small_traj.pdb"))
    topology = peng.to_openmm_topology(template)

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
        nonbondedCutoff=1.0*u.nanometers,
        constraints=app.HBonds,
        ewaldErrorTolerance=0.0005
    )


class testOpenMMEngine(object):
    def setUp(self):

        # OpenMM Integrator
        integrator = mm.LangevinIntegrator(
            300*u.kelvin,
            1.0/u.picoseconds,
            2.0*u.femtoseconds
        )
        integrator.setConstraintTolerance(0.00001)

        # Engine options
        options = {
            'nsteps_per_frame': 2,
            'platform': 'CPU',
            'solute_indices': range(22),
            'n_frames_max': 5,
            'timestep': 2.0 * u.femtoseconds
        }

        self.engine = peng.Engine(
            template,
            system,
            integrator,
            options
        )

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

        self.engine.current_snapshot = peng.Snapshot.construct(
            coordinates=np.array(testpos) * u.nanometers,
            box_vectors=np.zeros((3, 3)),
            velocities=np.array(testvel) * u.nanometers / u.picoseconds
        )
        state = self.engine.simulation.context.getState(getPositions=True,
                                                        getVelocities=True)
        sim_coords = state.getPositions(asNumpy=True) / u.nanometers
        sim_vels = state.getVelocities(asNumpy=True) / (u.nanometers/u.picoseconds)

        np.testing.assert_almost_equal(testpos, sim_coords, decimal=5)
        np.testing.assert_almost_equal(testvel, sim_vels, decimal=5)

    def test_generate_next_frame(self):
        snap0 = peng.Snapshot(
            statics=self.engine.current_snapshot.statics,
            kinetics=self.engine.current_snapshot.kinetics
        )
        new_snap = self.engine.generate_next_frame()
        assert(new_snap is not snap0)
        assert(new_snap.statics is not snap0.statics)
        assert(new_snap.kinetics is not snap0.kinetics)
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
        assert_equal(self.engine.snapshot_timestep, 4 * u.femtoseconds)
