"""
@author David W.H. Swenson
"""
from __future__ import division
from __future__ import absolute_import

import pytest
from builtins import range
from builtins import object
from past.utils import old_div
import numpy as np
from nose.tools import (assert_equal)
from nose.plugins.skip import SkipTest
try:
    import simtk.openmm as mm
    from simtk.openmm import app
except ImportError:
    mm = None
    app = None

import openpathsampling.engines.openmm as peng
import openpathsampling.engines as dyn

import openpathsampling as paths

from .test_helpers import (
    true_func, data_filename,
    assert_equal_array_array,
    assert_not_equal_array_array,
    raises_with_message_like, u, md,
    A2BEnsemble)

import logging
logging.getLogger('openpathsampling.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.ensemble').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.netcdfplus').setLevel(logging.CRITICAL)

def setup_module():
    global topology, template, system, nan_causing_template
    if not (u and mm and app and md):
        raise SkipTest
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

    # this is crude but does the trick
    nan_causing_template = template.copy()
    kinetics = template.kinetics.copy()
    # this is crude but does the trick
    kinetics.velocities = kinetics.velocities.copy()
    kinetics.velocities[0] = \
        (np.zeros(template.velocities.shape[1]) + 1000000.) * \
        u.nanometers / u.picoseconds
    nan_causing_template.kinetics = kinetics


class TestOpenMMEngine(object):
    def setup(self):

        # OpenMM Integrator
        integrator = mm.LangevinIntegrator(
            300*u.kelvin,
            old_div(1.0,u.picoseconds),
            2.0*u.femtoseconds
        )
        integrator.setConstraintTolerance(0.00001)

        # Engine options
        options = {
            'n_steps_per_frame': 2,
            'n_frames_max': 5
        }

        self.engine = peng.Engine(
            template.topology,
            system,
            integrator,
            options=options
        )

        self.engine.initialize('CPU')
        context = self.engine.simulation.context
        zero_array = np.zeros((template.topology.n_atoms, 3))
        context.setPositions(template.coordinates)
        context.setVelocities(u.Quantity(zero_array, old_div(u.nanometers, u.picoseconds)))

    def teardown(self):
        pass

    def test_sanity(self):
        assert_equal(self.engine.n_steps_per_frame, 2)
        assert_equal(self.engine.n_frames_max, 5)
        # TODO: add more sanity checks
        pass  # not quite a SkipTest, but a reminder to add more

    def test_snapshot_get(self):
        snap = self.engine.current_snapshot
        state = self.engine.simulation.context.getState(getVelocities=True,
                                                        getPositions=True)
        pos = old_div(state.getPositions(asNumpy=True), u.nanometers)
        vel = old_div(state.getVelocities(asNumpy=True), (old_div(u.nanometers, u.picoseconds)))
        assert_equal_array_array(old_div(snap.coordinates, u.nanometers), pos)
        assert_equal_array_array(old_div(snap.velocities, (old_div(u.nanometers, u.picoseconds))),
                                 vel)

    def test_snapshot_set(self):
        pdb_pos = (old_div(template.coordinates, u.nanometers))
        testvel = []
        testpos = []
        for i in range(len(pdb_pos)):
            testpos.append(list(np.array(pdb_pos[i]) + 
                                np.array([1.0, 1.0, 1.0]))
                          )
            testvel.append([0.1*i, 0.1*i, 0.1*i])

        testbvecs = 5.0 * np.identity(3)

        self.engine.current_snapshot = peng.Snapshot.construct(
            coordinates=np.array(testpos) * u.nanometers,
            box_vectors=np.array(testbvecs) * u.nanometers,
            velocities=np.array(testvel) * u.nanometers / u.picoseconds
        )
        state = self.engine.simulation.context.getState(getPositions=True,
                                                        getVelocities=True)
        sim_coords = old_div(state.getPositions(asNumpy=True), u.nanometers)
        sim_bvecs = old_div(state.getPeriodicBoxVectors(asNumpy=True), u.nanometers)
        sim_vels = old_div(state.getVelocities(asNumpy=True), (old_div(u.nanometers,u.picoseconds)))

        np.testing.assert_almost_equal(testpos, sim_coords, decimal=5)
        np.testing.assert_almost_equal(testbvecs, sim_bvecs, decimal=5)
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
        old_pos = old_div(snap0.coordinates, u.nanometers)
        new_pos = old_div(new_snap.coordinates, u.nanometers)
        old_vel = old_div(snap0.velocities, (old_div(u.nanometers, u.picoseconds)))
        new_vel = old_div(new_snap.velocities, (old_div(u.nanometers, u.picoseconds)))
        assert_equal(old_pos.shape, new_pos.shape)
        assert_equal(old_vel.shape, new_vel.shape)
        assert_not_equal_array_array(old_pos, new_pos)
        assert_not_equal_array_array(old_vel, new_vel)

    def test_generate(self):
        try:
            _ = self.engine.generate(self.engine.current_snapshot, [true_func])
        except dyn.EngineMaxLengthError as e:
            traj = e.last_trajectory
            assert_equal(len(traj), self.engine.n_frames_max)
        else:
            raise RuntimeError('Did not have correct MaxLengthError')

    def test_snapshot_timestep(self):
        assert_equal(self.engine.snapshot_timestep, 4 * u.femtoseconds)

    @raises_with_message_like(paths.engines.EngineMaxLengthError,
                              "Hit maximal length")
    def test_fail_length(self):
        self.engine.options['on_max_length'] = 'fail'
        self.engine.options['n_max_length'] = 2
        _ = self.engine.generate(self.engine.current_snapshot, [true_func])

    @raises_with_message_like(paths.engines.EngineMaxLengthError,
                              'Failed to generate trajectory without hitting '
                              'max length')
    def test_retry_length(self):
        self.engine.on_max_length = 'retry'
        self.engine.options['n_max_length'] = 2
        self.engine.options['retries_when_max_length'] = 2
        _ = self.engine.generate(self.engine.current_snapshot, [true_func])

    # OpenMM CPU will throw an error and not return a snapshot with nan
    @raises_with_message_like(dyn.EngineNaNError,
                              '`nan` in snapshot')
    def test_fail_nan(self):
        self.engine.on_max_length = 'retry'
        self.engine.options['n_max_length'] = 2
        self.engine.options['retries_when_max_length'] = 2
        _ = self.engine.generate(nan_causing_template, [true_func])

    def test_nan_rejected(self):
        stateA = paths.EmptyVolume()  # will run indefinitely
        stateB = paths.EmptyVolume()
        tps = A2BEnsemble(stateA, stateB)
        self.engine.n_frames_max = 10

        init_traj = paths.Trajectory([nan_causing_template] * 5)
        init_samp = paths.SampleSet([paths.Sample(
            trajectory=init_traj,
            replica=0,
            ensemble=tps
        )])

        mover = paths.BackwardShootMover(
            ensemble=tps,
            selector=paths.UniformSelector(),
            engine=self.engine
        )
        change = mover.move(init_samp)

        assert (isinstance(change, paths.RejectedNaNSampleMoveChange))
        assert_equal(change.details.rejection_reason, 'nan')
        # since we shoot, we start with a shorter trajectory
        assert(len(change.samples[0].trajectory) < len(init_traj))

        newsamp = init_samp.apply_samples(change)
        assert_equal(len(newsamp), 1)

        # make sure there is no change!
        assert_equal(init_samp[0].trajectory, init_traj)

    def test_max_length_rejected(self):
        stateA = paths.EmptyVolume()  # will run indefinitely
        stateB = paths.EmptyVolume()
        tps = A2BEnsemble(stateA, stateB)
        self.engine.options['n_frames_max'] = 10
        self.engine.on_max_length = 'fail'

        init_traj = paths.Trajectory([template] * 5)
        init_samp = paths.SampleSet([paths.Sample(
            trajectory=init_traj,
            replica=0,
            ensemble=tps
        )])

        mover = paths.BackwardShootMover(
            ensemble=tps,
            selector=paths.UniformSelector(),
            engine=self.engine
        )
        change = mover.move(init_samp)

        assert(isinstance(change, paths.RejectedMaxLengthSampleMoveChange))
        assert_equal(change.details.rejection_reason, 'max_length')
        assert_equal(
            len(change.samples[0].trajectory), self.engine.n_frames_max)

        newsamp = init_samp.apply_samples(change)
        assert_equal(len(newsamp), 1)

        # make sure there is no change!
        assert_equal(init_samp[0].trajectory, init_traj)

    @pytest.mark.parametrize('has_constraints', [True, False])
    def test_has_constraints(self, has_constraints):
        omt = pytest.importorskip('openmmtools')
        constraints = app.HBonds if has_constraints else None
        system = omt.testsystems.AlanineDipeptideVacuum(
            constraints=constraints
        )
        template = peng.snapshot_from_testsystem(system)
        engine = peng.Engine(
            topology=template.topology,
            system=system.system,
            integrator=omt.integrators.VVVRIntegrator()
        )
        assert engine.has_constraints() == has_constraints
