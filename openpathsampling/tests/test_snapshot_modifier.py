from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        raises, assert_almost_equal, assert_true)
from nose.plugins.skip import SkipTest
from numpy.testing import assert_array_almost_equal
from test_helpers import make_1d_traj

import openpathsampling as paths
import openpathsampling.engines as peng
import numpy as np

import openmmtools as omt
import simtk.unit as u
import openpathsampling.engines.openmm as omm_engine

from openpathsampling.snapshot_modifier import *

import logging
logging.getLogger('openpathsampling.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.netcdfplus').setLevel(logging.CRITICAL)

class testNoModification(object):
    def setup(self):
        self.modifier = NoModification()
        self.snapshot_1D = peng.toy.Snapshot(
            coordinates=np.array([0.0, 1.0, 2.0, 3.0]),
            velocities=np.array([0.5, 1.5, 2.5, 3.5])
        )
        self.snapshot_3D = peng.openmm.MDSnapshot(
            coordinates=np.array([[0.0, 0.1, 0.2],
                                  [1.0, 1.1, 1.2],
                                  [2.0, 2.1, 2.2],
                                  [3.0, 3.1, 3.2]]),
            velocities=np.array([[0.5, 0.6, 0.7],
                                 [1.5, 1.6, 1.7],
                                 [2.5, 2.6, 2.7],
                                 [3.5, 3.6, 3.7]])
        )

    # test methods from the abstract base class along with NoModification
    def test_extract_subset(self):
        mod = NoModification(subset_mask=[1,2])
        sub_1Dx = mod.extract_subset(self.snapshot_1D.coordinates)
        assert_array_almost_equal(sub_1Dx, np.array([1.0, 2.0]))
        sub_1Dv = mod.extract_subset(self.snapshot_1D.velocities)
        assert_array_almost_equal(sub_1Dv, np.array([1.5, 2.5]))

        sub_3Dx = mod.extract_subset(self.snapshot_3D.coordinates)
        assert_array_almost_equal(sub_3Dx, np.array([[1.0, 1.1, 1.2],
                                                     [2.0, 2.1, 2.2]]))
        sub_3Dv = mod.extract_subset(self.snapshot_3D.velocities)
        assert_array_almost_equal(sub_3Dv, np.array([[1.5, 1.6, 1.7],
                                                     [2.5, 2.6, 2.7]]))

    def test_apply_to_subset(self):
        mod = NoModification(subset_mask=[1,2])
        copy_1Dx = self.snapshot_1D.coordinates.copy()
        new_1Dx = mod.apply_to_subset(copy_1Dx, np.array([-1.0, -2.0]))
        assert_array_almost_equal(new_1Dx, np.array([0.0, -1.0, -2.0, 3.0]))
        # and check that memory points to the right things; orig unchanged
        assert_true(copy_1Dx is new_1Dx)
        assert_array_almost_equal(self.snapshot_1D.coordinates,
                                  np.array([0.0, 1.0, 2.0, 3.0]))

        copy_3Dx = self.snapshot_3D.coordinates.copy()
        new_3Dx = mod.apply_to_subset(copy_3Dx,
                                      np.array([[-1.0, -1.1, -1.2],
                                                [-2.0, -2.1, -2.2]]))
        assert_array_almost_equal(new_3Dx, np.array([[0.0, 0.1, 0.2],
                                                     [-1.0, -1.1, -1.2],
                                                     [-2.0, -2.1, -2.2],
                                                     [3.0, 3.1, 3.2]]))
        # and check that memory points to the right things; orig unchanged
        assert_true(copy_3Dx is new_3Dx)
        assert_array_almost_equal(self.snapshot_3D.coordinates,
                                  np.array([[0.0, 0.1, 0.2],
                                            [1.0, 1.1, 1.2],
                                            [2.0, 2.1, 2.2],
                                            [3.0, 3.1, 3.2]]))

    def test_call(self):
        new_1D = self.modifier(self.snapshot_1D)
        assert_array_almost_equal(self.snapshot_1D.coordinates,
                                  new_1D.coordinates)
        assert_array_almost_equal(self.snapshot_1D.velocities,
                                  new_1D.velocities)
        new_3D = self.modifier(self.snapshot_3D)
        assert_array_almost_equal(self.snapshot_3D.coordinates,
                                  new_3D.coordinates)
        assert_array_almost_equal(self.snapshot_3D.velocities,
                                  new_3D.velocities)
        assert_true(self.snapshot_1D.coordinates is not new_1D.coordinates)
        assert_true(self.snapshot_1D.velocities is not new_1D.velocities)
        assert_true(self.snapshot_3D.coordinates is not new_3D.coordinates)
        assert_true(self.snapshot_3D.velocities is not new_3D.velocities)


class testRandomizeVelocities(object):
    def setup(self):
        # TODO: check against several possibilities, including various
        # combinations of shapes of velocities and masses.
        topology_2x3D = paths.engines.toy.Topology(
            n_spatial=3, n_atoms=2, masses=np.array([2.0, 3.0]), pes=None
        )
        topology_3x1D = paths.engines.toy.Topology(
            n_spatial=1, n_atoms=3, masses=np.array([[2.0], [3.0], [4.0]]),
            pes=None
        )
        topology_1x2D = paths.engines.toy.Topology(
            n_spatial=2, n_atoms=1, masses=np.array([1.0, 2.0]), pes=None
        )
        self.snap_2x3D = paths.engines.toy.Snapshot(
            coordinates=np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]),
            velocities=np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]),
            engine=paths.engines.toy.Engine({}, topology_2x3D)
        )
        self.snap_3x1D = paths.engines.toy.Snapshot(
            coordinates=np.array([[0.0], [0.0], [0.0]]),
            velocities=np.array([[0.0], [0.0], [0.0]]),
            engine=paths.engines.toy.Engine({}, topology_3x1D)
        )
        self.snap_1x2D = paths.engines.toy.Snapshot(
            coordinates=np.array([[0.0, 0.0]]),
            velocities=np.array([[0.0, 0.0]]),
            engine=paths.engines.toy.Engine({}, topology_1x2D)
        )

    def test_call(self):
        # NOTE: these tests basically check the API. Tests for correctness
        # are in `test_snapshot_modifier.ipynb`, because they are inherently
        # stochastic.
        randomizer = RandomVelocities(beta=1.0/5.0)
        new_1x2D = randomizer(self.snap_1x2D)
        assert_equal(new_1x2D.coordinates.shape, new_1x2D.velocities.shape)
        assert_array_almost_equal(new_1x2D.coordinates,
                                  self.snap_1x2D.coordinates)
        assert_true(new_1x2D is not self.snap_1x2D)
        assert_true(new_1x2D.coordinates is not self.snap_1x2D.coordinates)
        assert_true(new_1x2D.velocities is not self.snap_1x2D.velocities)
        for val in new_1x2D.velocities.flatten():
            assert_not_equal(val, 0.0)

        new_2x3D = randomizer(self.snap_2x3D)
        assert_equal(new_2x3D.coordinates.shape, new_2x3D.velocities.shape)
        assert_array_almost_equal(new_2x3D.coordinates,
                                  self.snap_2x3D.coordinates)
        assert_true(new_2x3D is not self.snap_2x3D)
        assert_true(new_2x3D.coordinates is not self.snap_2x3D.coordinates)
        assert_true(new_2x3D.velocities is not self.snap_2x3D.velocities)
        for val in new_2x3D.velocities.flatten():
            assert_not_equal(val, 0.0)

        new_3x1D = randomizer(self.snap_3x1D)
        assert_equal(new_3x1D.coordinates.shape, new_3x1D.velocities.shape)
        assert_array_almost_equal(new_3x1D.coordinates,
                                  self.snap_3x1D.coordinates)
        assert_true(new_3x1D is not self.snap_3x1D)
        assert_true(new_3x1D.coordinates is not self.snap_3x1D.coordinates)
        assert_true(new_3x1D.velocities is not self.snap_3x1D.velocities)
        for val in new_3x1D.velocities.flatten():
            assert_not_equal(val, 0.0)

    def test_subset_call(self):
        randomizer = RandomVelocities(beta=1.0/5.0, subset_mask=[0])
        new_2x3D = randomizer(self.snap_2x3D)
        assert_equal(new_2x3D.coordinates.shape, new_2x3D.velocities.shape)
        assert_array_almost_equal(new_2x3D.coordinates,
                                  self.snap_2x3D.coordinates)
        assert_true(new_2x3D is not self.snap_2x3D)
        assert_true(new_2x3D.coordinates is not self.snap_2x3D.coordinates)
        assert_true(new_2x3D.velocities is not self.snap_2x3D.velocities)
        # show that the unchanged atom is, in fact, unchanged
        assert_array_almost_equal(new_2x3D.velocities[1],
                                  self.snap_2x3D.velocities[1])
        for val in new_2x3D.velocities[0]:
            assert_not_equal(val, 0.0)

    def test_with_openmm_snapshot(self):
        # note: this is only a smoke test; correctness depends on OpenMM's
        # tests of its constraint approaches.
        test_system = omt.testsystems.AlanineDipeptideVacuum()
        template = omm_engine.snapshot_from_testsystem(test_system)
        engine = omm_engine.Engine(
            topology=template.topology,
            system=test_system.system,
            integrator=omt.integrators.VVVRIntegrator()
        )
        beta = 1.0 / (300.0 * u.kelvin * u.BOLTZMANN_CONSTANT_kB)

        # when the engine doesn't have an existing snapshot
        randomizer = RandomVelocities(beta=beta, engine=engine)
        new_snap = randomizer(template)
        # coordinates stayed the same
        assert_array_almost_equal(template.coordinates,
                                  new_snap.coordinates)
        # velocities changed
        assert_equal(np.isclose(template.velocities,
                                new_snap.velocities).all(),
                     False)
        engine.generate(new_snap, [lambda x, foo: len(x) <= 4])

        # when the engine does have an existing snapshot
        zeros = np.zeros((engine.n_atoms, engine.n_spatial))
        zero_snap = paths.engines.openmm.Snapshot.construct(
            coordinates=zeros * u.nanometer,
            velocities=zeros * u.nanometer / u.picosecond,
            box_vectors=template.box_vectors,
            engine=engine
        )
        engine.current_snapshot = zero_snap
        randomizer = RandomVelocities(beta=beta, engine=engine)
        new_snap = randomizer(template)
        # coordinates stayed the same
        assert_array_almost_equal(template.coordinates,
                                  new_snap.coordinates)
        # velocities changed
        assert_equal(np.isclose(template.velocities,
                                new_snap.velocities).all(),
                     False)
        # internal snapshot unchanged
        assert_equal(engine.current_snapshot, zero_snap)
        engine.generate(new_snap, [lambda x, foo: len(x) <= 4])

class testGeneralizedDirectionModifier(object):
    def setup(self):
        import openpathsampling.engines.toy as toys
        # applies one delta_v to all atoms
        self.toy_modifier_all = GeneralizedDirectionModifier(1.5)
        # defines delta_v per atom, including those not in the mask
        self.toy_modifier_long_dv = GeneralizedDirectionModifier(
            delta_v=[0.5, 1.0, 2.0],
            subset_mask=[1, 2]
        )
        # defines delta_v per atom in the subset mask
        self.toy_modifier = GeneralizedDirectionModifier(
            delta_v=[1.0, 2.0],
            subset_mask=[1, 2]
        )
        self.toy_engine = toys.Engine(
            topology=toys.Topology(n_spatial=2, n_atoms=3, pes=None,
                                   masses=[1.0, 1.5, 4.0]),
            options={}
        )
        self.toy_snapshot = toys.Snapshot(
            coordinates=np.array([[0.0, 0.0], [0.0, 0.0], [0.0, 0.0]]),
            velocities=np.array([[1.0, 1.0], [2.0, 2.0], [3.0, 3.0]]),
            engine=self.toy_engine
        )

        # create the OpenMM versions
        self.openmm_modifier = GeneralizedDirectionModifier(1.0)
        ad_vacuum = omt.testsystems.AlanineDipeptideVacuum(constraints=None)
        self.test_snap = omm_engine.snapshot_from_testsystem(ad_vacuum)
        self.openmm_engine = omm_engine.Engine(
            topology=self.test_snap.topology,
            system=ad_vacuum.system,
            integrator=omt.integrators.VVVRIntegrator()
        )
        self.openmm_snap = self.test_snap.copy_with_replacement(
            engine=self.openmm_engine
        )
        pass

    def test_verify_snapshot_toy(self):
        self.toy_modifier._verify_snapshot(self.toy_snapshot)
        self.toy_modifier_all._verify_snapshot(self.toy_snapshot)
        self.toy_modifier_long_dv._verify_snapshot(self.toy_snapshot)

    def test_verify_snapshot_openmm(self):
        self.openmm_modifier._verify_snapshot(self.openmm_snap)

    @raises(RuntimeError)
    def test_verify_snapshot_no_dofs(self):
        assert_true(isinstance(self.test_snap.engine,
                               omm_engine.tools.OpenMMToolsTestsystemEngine))
        self.openmm_modifier._verify_snapshot(self.test_snap)

    @raises(RuntimeError)
    def test_verify_snapshot_constraints(self):
        ad_vacuum_constr = omt.testsystems.AlanineDipeptideVacuum()
        constrained_engine = omm_engine.Engine(
            topology=self.test_snap.topology,
            system=ad_vacuum_constr.system,
            integrator=omt.integrators.VVVRIntegrator()
        )
        constr_snap = self.test_snap.copy_with_replacement(
            engine=constrained_engine
        )
        self.openmm_modifier._verify_snapshot(constr_snap)

    def test_verify_snapshot_box_vectors(self):
        ad_explicit = omt.testsystems.AlanineDipeptideExplicit(
            constraints=None,
            rigid_water=False
        )
        ad_explicit_tmpl = omm_engine.snapshot_from_testsystem(ad_explicit)
        explicit_engine= omm_engine.Engine(
            topology=ad_explicit_tmpl.topology,
            system=ad_explicit.system,
            integrator=omt.integrators.VVVRIntegrator()
        )
        ad_explicit_snap = ad_explicit_tmpl.copy_with_replacement(
            engine=explicit_engine
        )
        self.openmm_modifier._verify_snapshot(ad_explicit_snap)

    def test_dv_widths_toy(self):
        pass

    def test_dv_widths_openmm(self):
        pass
