'''
@author: David W.H. Swenson
'''

import os
import numpy as np

from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        assert_almost_equal, raises)
from nose.plugins.skip import Skip, SkipTest
from test_helpers import true_func, assert_equal_array_array

from openpathsampling.toy_dynamics.toy_pes import *
from openpathsampling.toy_dynamics.toy_integrators import *
from openpathsampling.toy_dynamics.toy_engine import *
from openpathsampling.topology import ToyTopology
from openpathsampling.snapshot import Snapshot, Momentum, Configuration


# =========================================================================
# This single test module includes all the tests for the toy_dynamics
# subpackage. 
# =========================================================================

def setUp():
    # set up globals
    global gaussian, linear, outer, harmonic
    gaussian = Gaussian(6.0, [2.5, 40.0], [0.8, 0.5])
    outer = OuterWalls([1.12, 2.0], [0.2, -0.25])
    linear = LinearSlope([1.5, 0.75], 0.5)
    harmonic = HarmonicOscillator([1.5, 2.0], [0.5, 3.0], [0.25, 0.75])
    global init_pos, init_vel, sys_mass
    init_pos = np.array([0.7, 0.65])
    init_vel = np.array([0.6, 0.5])
    sys_mass = np.array([1.5, 1.5])

# === TESTS FOR TOY POTENTIAL ENERGY SURFACES =============================

class testHarmonicOscillator(object):
    def setUp(self):
        self.positions = init_pos
        self.velocities = init_vel
        self.mass = sys_mass

    def test_V(self):
        # k = m * omega^2 = [1.5, 1.5] * [0.5, 3.0]^2
        #                 = [0.375, 13.5]
        # V = 0.5*( 1.5*0.375*((0.7)-0.25)^2 + 2.0*13.5*((0.65)-0.75)^2)
        #   = 0.191953125
        assert_almost_equal(harmonic.V(self), 0.191953125)

    def test_dVdx(self):
        # [1.5, 2.0] * [1.5, 1.5] * [0.5, 3.0]^2 * [(0.7)-0.25, (0.65)-0.75]
        #   = [1.5*1.5*0.5^2*((0.7)-0.25), 2.0*1.5*3.0^2*((0.65)-0.75)]
        #   = [0.253125, -2.7]
        for (experiment, theory) in zip(harmonic.dVdx(self),
                                        [0.253125, -2.7]):
            assert_almost_equal(experiment, theory)

class testGaussian(object):
    def setUp(self):
        self.positions = init_pos
        self.velocities = init_vel
        self.mass = sys_mass

    def test_V(self):
        # 6.0*exp(-2.5*((0.7)-0.8)^2-40.0*((0.65)-0.5)^2) = 2.37918851445
        assert_almost_equal(gaussian.V(self), 2.37918851445)

    def test_dVdx(self):
        # exp(-2.5*((0.7)-0.8)^2-40*((0.65)-0.5)^2)*(-30*(0.7)+24)
        assert_almost_equal(gaussian.dVdx(self)[0], 1.18959425722)
        # -480*((0.65)-0.5)*exp(-2.5*((0.7)-0.8)^2-40*((0.65)-0.5)^2)
        assert_almost_equal(gaussian.dVdx(self)[1], -28.5502621734)

class testOuterWalls(object):
    def setUp(self):
        self.positions = init_pos
        self.velocities = init_vel
        self.mass = sys_mass

    def test_V(self):
        # 1.12*(0.7-0.2)^6+2.0*(0.65-(-0.25))^6 = 1.080382
        assert_almost_equal(outer.V(self), 1.080382)

    def test_dVdx(self):
        # 6*1.12*(0.7-0.2)^5 = 0.21
        assert_almost_equal(outer.dVdx(self)[0], 0.21)
        # 6*2.0*(0.65-(-0.25))^5 = 7.08588
        assert_almost_equal(outer.dVdx(self)[1], 7.08588)

class testLinearSlope(object):
    def setUp(self):
        self.positions = init_pos
        self.velocities = init_vel
        self.mass = sys_mass

    def test_V(self):
        assert_almost_equal(linear.V(self), 2.0375)

    def test_dVdx(self):
        assert_equal(linear.dVdx(self), [1.5, 0.75])

class testCombinations(object):
    def setUp(self):
        self.positions = init_pos
        self.velocities = init_vel
        self.mass = sys_mass
        self.simpletest = gaussian + gaussian
        self.fullertest = gaussian + outer - linear

    def test_V(self):
        assert_almost_equal(self.simpletest.V(self), 2*2.37918851445)
        assert_almost_equal(self.fullertest.V(self), 
                            2.37918851445 + 1.080382 - 2.0375)

    def test_dVdx(self):
        assert_almost_equal(self.simpletest.dVdx(self)[0], 2*1.18959425722)
        assert_almost_equal(self.simpletest.dVdx(self)[1], 2*-28.5502621734)
        assert_almost_equal(self.fullertest.dVdx(self)[0],
                            1.18959425722 + 0.21 - 1.5)
        assert_almost_equal(self.fullertest.dVdx(self)[1],
                            -28.5502621734 + 7.08588 - 0.75)

    def test_kinetic_energy(self):
        assert_almost_equal(self.simpletest.kinetic_energy(self), 0.4575)


# === TESTS FOR TOY ENGINE OBJECT =========================================

class test_convert_fcn(object):
    def test_convert_to_3Ndim(v):
        assert_equal_array_array(convert_to_3Ndim([1.0, 2.0]),
                                 np.array([[1.0, 2.0, 0.0]]))
        assert_equal_array_array(convert_to_3Ndim([1.0, 2.0, 3.0]), 
                                 np.array([[1.0, 2.0, 3.0]]))
        assert_equal_array_array(convert_to_3Ndim([1.0, 2.0, 3.0, 4.0]),
                                 np.array([[1.0, 2.0, 3.0], [4.0, 0.0, 0.0]]))

class testToyEngine(object):
    def setUp(self):
        pes = linear
        integ = LeapfrogVerletIntegrator(dt=0.002)
        topology=ToyTopology(
            n_spatial = 2,
            masses = sys_mass,
            pes = pes
        )
        template = Snapshot(
            coordinates=init_pos.copy(),
            velocities=init_pos.copy(),
            potential_energy = 0.0,
            kinetic_energy = 0.0,
            topology=topology
        )
        options={
            'integ' : integ,
            'n_frames_max' : 5}
        sim = ToyEngine(options=options,
                        template=template
                       )

        sim.positions = init_pos.copy()
        sim.velocities = init_vel.copy()
        sim.nsteps_per_frame = 10
        self.sim = sim

    def teardown(self):
        if os.path.isfile('toy_tmp.nc'):
            os.remove('toy_tmp.nc')

    def test_sanity(self):
        assert_items_equal(self.sim.mass, sys_mass)
        assert_items_equal(self.sim._minv, [1.0/m_i for m_i in sys_mass])
        assert_equal(self.sim.nsteps_per_frame, 10)

    def test_momentum_getter(self):
        momentum = self.sim.momentum
        assert_items_equal(momentum.velocities[0],
                           self.sim.velocities)
        assert_equal(momentum.kinetic_energy,
                     self.sim.pes.kinetic_energy(self.sim))
    
    def test_momentum_setter(self):
        raise SkipTest()
        self.sim.momentum = Momentum(velocities=np.array([[4, 5, 6]]))
        assert_items_equal(self.sim.velocities, [4, 5])

    def test_configuration_setter(self):
        raise SkipTest()
        # This should not work anymore
        self.sim.configuration = Configuration(
            coordinates=np.array([[1, 2, 3]])
        )
        assert_items_equal(self.sim.positions, [1, 2])


    def test_load_configuration(self):
        configuration = self.sim.configuration
        assert_items_equal(configuration.coordinates[0],
                           self.sim.positions)
        assert_equal(configuration.potential_energy,
                     self.sim.pes.V(self.sim))
        assert_equal(configuration.box_vectors, None)

    def test_snapshot_get(self):
        snapshot = self.sim.current_snapshot
        n_spatial=self.sim.n_spatial
        assert_items_equal(snapshot.momentum.velocities[0],
                           self.sim.velocities)
        assert_equal(snapshot.momentum.kinetic_energy,
                     self.sim.pes.kinetic_energy(self.sim))
        assert_items_equal(snapshot.configuration.coordinates[0],
                           self.sim.positions)
        assert_equal(snapshot.configuration.potential_energy,
                     self.sim.pes.V(self.sim))
        assert_equal(snapshot.configuration.box_vectors, None)
        
    def test_snapshot_set(self):
        raise SkipTest()
        snap = Snapshot(coordinates=np.array([[1,2,3]]), 
                        velocities=np.array([[4,5,6]]))
        # note that we truncate the z-direction, regardless of what it is.
        # This is for a 2D model!
        self.sim.current_snapshot = snap
        assert_items_equal(self.sim.positions, [1,2])
        assert_items_equal(self.sim.velocities, [4,5])

    def test_generate_next_frame(self):
        # we test correctness by integrating forward, then backward
        assert_items_equal(self.sim.positions, init_pos)
        assert_items_equal(self.sim.velocities, init_vel)
        snap = self.sim.generate_next_frame()
        #assert_equal_array_array(snap.coordinates,
                                 #np.array([init_pos.append(0.0)]))
        self.sim.velocities = -self.sim.velocities
        snap2 = self.sim.generate_next_frame()
        np.testing.assert_allclose(snap2.coordinates[0], init_pos)

    def test_generate(self):
        self.sim.initialized = True
        traj = self.sim.generate(self.sim.current_snapshot, [true_func])
        assert_equal(len(traj), self.sim.n_frames_max)

    @raises(RuntimeWarning)
    def test_generate_uninitialized(self):
        traj = self.sim.generate(self.sim.current_snapshot, [true_func])

    def test_start_with_snapshot(self):
        snap = Snapshot(coordinates=np.array([1,2]), 
                        velocities=np.array([3,4]))
        self.sim.start(snapshot=snap)
        self.sim.stop([snap])



# === TESTS FOR TOY INTEGRATORS ===========================================

class testLeapfrogVerletIntegrator(object):
    def setUp(self):
        pes = linear
        integ = LeapfrogVerletIntegrator(dt=0.002)
        topology=ToyTopology(
            n_spatial = 2,
            masses = sys_mass,
            pes = pes
        )
        template = Snapshot(
            coordinates=init_pos.copy(),
            velocities=init_pos.copy(),
            potential_energy = 0.0,
            kinetic_energy = 0.0,
            topology=topology
        )
        options={
            'integ' : integ,
            'n_frames_max' : 5}
        sim = ToyEngine(options=options,
                        template=template
                       )

        sim.positions = init_pos.copy()
        sim.velocities = init_vel.copy()

        sim.nsteps_per_frame = 10
        self.sim = sim

    def test_momentum_update(self):
        self.sim.integ._momentum_update(self.sim, 0.002)
        # velocities = init_vel - pes.dVdx(init_pos)/m*dt 
        #            = [0.6, 0.5] - [1.5, 0.75]/[1.5, 1.5] * 0.002
        #            = [0.598, 0.499]
        assert_equal(self.sim.velocities[0], 0.598)
        assert_equal(self.sim.velocities[1], 0.499)
    
    def test_position_update(self):
        self.sim.integ._position_update(self.sim, 0.002)
        # positions = init_pos + velocities * dt
        #           = [0.7, 0.65] + [0.6, 0.5]*0.002 = [0.7012, 0.651]
        assert_almost_equal(self.sim.positions[0], 0.7012)
        assert_almost_equal(self.sim.positions[1], 0.651)

    def test_step(self):
        # no assertions since the tests of position/momentum updates should
        # handle that... this is just to make sure we run
        self.sim.integ.step(self.sim, 2)



class testLangevinBAOABIntegrator(object):
    '''Testing for correctness is hard, since this is a stochastic
    calculation. However, we can at least run tests to make sure nothing
    crashes.'''
    def setUp(self):
        pes = linear
        integ = LangevinBAOABIntegrator(dt=0.002, temperature=0.5,
                                        gamma=1.0)
        topology=ToyTopology(
            n_spatial = 2,
            masses = sys_mass,
            pes = pes
        )
        template = Snapshot(
            coordinates=init_pos.copy(),
            velocities=init_pos.copy(),
            potential_energy = 0.0,
            kinetic_energy = 0.0,
            topology=topology
        )
        options={
            'integ' : integ,
            'n_frames_max' : 5}
        sim = ToyEngine(options=options,
                        template=template
                       )

        sim.positions = init_pos.copy()
        sim.velocities = init_vel.copy()

        sim.nsteps_per_frame = 10
        self.sim = sim

    def test_OU_update(self):
        # we can't actually test for correct values, but we *can* test that
        # some things *don't* happen
        assert_equal(self.sim.velocities[0], init_vel[0])
        assert_equal(self.sim.velocities[1], init_vel[1])
        self.sim.integ._OU_update(self.sim, 0.002)
        assert_not_equal(self.sim.velocities[0], init_vel[0])
        assert_not_equal(self.sim.velocities[1], init_vel[1])
        # tests that the same random number wasn't used for both:
        assert_not_equal(self.sim.velocities[0] - init_vel[0],
                         self.sim.velocities[1] - init_vel[1])

    def test_step(self):
        self.sim.generate_next_frame()
