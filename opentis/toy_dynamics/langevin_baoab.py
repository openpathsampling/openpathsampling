import math

class LangevinBAOAB(object):
    """ Langevin integrator for simple toy models

    Implementation of the BAOAB integrator by Matthews and Leimkuhler.

    Reference
    ---------

    """

    def __init__(self, dt, nt, temperature, gamma):
        self.beta = 1.0 / temperature
        self.c1 = math.exp(-gamma*dt)
        self.c2 = (1.0-self.c1)/gamma
        self.c3 = math.sqrt((1.0 - self.c1*self.c1) / self.beta)

    def _momentum_update(self, config):
        pdot = self.pes.pdot(config)
        pass

    def _position_update(self, config):
        pass

    def _OU_update(self, config):
        pass

    def step(self, config):
        pass
