import numpy as np
from nose.tools import (assert_almost_equal)
from openpathsampling.engines.toy.sshooting_integrators import OverdampedLangevinIntegrator
import openpathsampling.engines.toy as toys

def setUp():
    global linear
    linear = toys.LinearSlope([1.5, 0.75], 0.5)
    global init_pos, init_vel, sys_mass
    init_pos = np.array([0.7, 0.65])
    init_vel = np.array([0.6, 0.5])
    sys_mass = np.array([1.5, 1.5])

class testOverdampedLangevinIntegrator(object):
    '''This is only a test if the integrator runs. Because of its stochastic
    nature we can not actually test its correctness easily.'''
    def setUp(self):
        pes = linear
        integ = OverdampedLangevinIntegrator(dt=0.001, temperature=4.0, D=1.0)
        topology=toys.Topology(
            n_spatial = 2,
            masses = sys_mass,
            pes = pes
        )
        options={
            'integ' : integ,
            'n_frames_max' : 5}
        sim = toys.Engine(options=options, topology=topology)

        template = toys.Snapshot(
            coordinates=init_pos.copy(),
            velocities=init_vel.copy(),
            engine=sim
        )

        sim.positions = init_pos.copy()
        sim.velocities = init_vel.copy()

        sim.n_steps_per_frame = 10
        self.sim = sim

    def test_step(self):
        self.sim.generate_next_frame()
