'''
@author: David W.H. Swenson
'''
from nose.tools import assert_equal, assert_almost_equal, raises
from nose.plugins.skip import Skip, SkipTest

from opentis.toy_dynamics.toy_pes import *
from opentis.toy_dynamics.toy_integrators import *
from opentis.toy_dynamics.toy_simulation import *
#from opentis.snapshot import Snapshot

# =========================================================================
# This single test module includes all the tests for the toy_dynamics
# subpackage. 
# =========================================================================

def setUp():
    # set up globals
    global gaussian, linear, outer
    gaussian = Gaussian(6.0, [2.5, 40.0], [0.8, 0.5])
    outer = OuterWalls([1.12, 2.0], [0.2, -0.25])
    linear = LinearSlope([1.5, 0.75], 0.5)
    global init_pos, init_vel, sys_mass
    init_pos = [0.7, 0.65]
    init_vel = [0.6, 0.5]
    sys_mass = [1.5, 1.5]

# === TESTS FOR TOY POTENTIAL ENERGY SURFACES =============================

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


# === TESTS FOR TOY SIMULATION OBJECT =====================================

class testToySimulation(object):
    def setUp(self):
        pes = linear
        integ = EulerIntegrator(dt=0.002)
        sim = ToySimulation(pes, integ)
        sim.positions = init_pos
        sim.velocities = init_vel
        sim.mass = sys_mass
        sim.nsteps_per_iteration = 10
        self.sim = sim

    def test_sanity(self):
        assert_equal(self.sim.mass, sys_mass)
        assert_equal(self.sim.minv[0], 1.0/sys_mass[0])
        assert_equal(self.sim.minv[1], 1.0/sys_mass[1])
        assert_equal(self.sim.nsteps_per_iteration, 10)

    def test_load_momentum(self):
        raise SkipTest

    def test_load_configuration(self):
        raise SkipTest

    def test_load_snapshot(self):
        raise SkipTest

    def test_init_simulation_with_snapshot(self):
        raise SkipTest

    def test_generate_next_frame(self):
        raise SkipTest



# === TESTS FOR TOY INTEGRATORS ===========================================

class testEulerIntegrator(object):
    def setUp(self):
        pes = linear
        integ = EulerIntegrator(dt=0.002)
        sim = ToySimulation(pes, integ)
        sim.positions = init_pos
        sim.velocities = init_vel
        sim.mass = sys_mass
        sim.nsteps_per_iteration = 10
        self.sim = sim
    
    def test_momentum_update(self):
        raise SkipTest
    
    def test_position_update(self):
        raise SkipTest

    def test_step(self):
        raise SkipTest

class testLangevinBAOABIntegrator(object):
    '''Testing for correctness is hard, since this is a stochastic
    calculation. However, we can at least run tests to make sure nothing
    crashes.'''
    def setUp(self):
        pes = linear
        integ = LangevinBAOABIntegrator(dt=0.002, temperature=0.5,
                                        gamma=1.0)
        sim = ToySimulation(pes, integ)
        sim.positions = init_pos
        sim.velocities = init_vel
        sim.mass = sys_mass
        sim.nsteps_per_iteration = 10
        self.sim = sim

    def test_step(self):
        self.sim.generate_next_frame()
