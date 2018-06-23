import numpy as np
from nose.tools import (assert_almost_equal)
from openpathsampling.engines.toy.sshooting_pes import DoubleWell

def setUp():
    # set up globals
    global doublewell
    doublewell = DoubleWell([5.0, 2.0], [2.0, 3.0])

    global init_pos, init_vel, sys_mass
    init_pos = np.array([0.7, 0.65])
    init_vel = np.array([0.6, 0.5])
    sys_mass = np.array([1.5, 1.5])

class testDoubleWell(object):
    def setUp(self):
        self.positions = init_pos
        self.velocities = init_vel
        self.mass = sys_mass

    def test_V(self):
        assert_almost_equal(doublewell.V(self), 208.7475125)

    def test_dVdx(self):
        for (experiment, theory) in zip(doublewell.dVdx(self),
                                        [-49.14, -44.603]):
            assert_almost_equal(experiment, theory)
