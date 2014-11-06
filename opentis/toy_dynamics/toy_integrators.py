import math

from numpy.random import normal as random_normal

class EulerIntegrator(object):
    """ Euler Integrator. Not for actual use; but the momentum and position
    update functions are used in other integrators, so we inherit from this.
    """
    def __init__(self, dt):
        self.dt = dt

    def _momentum_update(self, sys, mydt):
        sys.velocities -= sys.pes.dVdx(sys) * sys.minv * mydt

    def _position_update(self, sys, mydt):
        sys.positions += sys.pes.velocities * mydt

    def step(self, sys):
        self._momentum_update(sys, mydt)
        self._position_update(sys, mydt)


class LangevinBAOABIntegrator(EulerIntegrator):
    """ Langevin integrator for simple toy models

    Implementation of the BAOAB integrator of Leimkuhler and Matthews. In
    particular, see the appendix on p.54 of the reference below, which is
    where we take our notation from.

    Reference
    ---------
    B. Leimkuhler and C. Matthews. "Rational Construction of Stochastic
    Numerical Methods for Molecular Sampling." Appl. Math. Res. Express,
    2013, 34-56 (2013). doi:10.1093/amrx/abs010

    """

    def __init__(self, dt, temperature, gamma):
        self.beta = 1.0 / temperature
        self.c1 = math.exp(-gamma*dt)
        self.c2 = (1.0-self.c1)/gamma
        self.c3 = math.sqrt((1.0 - self.c1*self.c1) / self.beta)
        self.dt = dt

    def _OU_update(self, sys, mydt):
        R = random_normal(size=len(sys.velocities))
        self.velocities = (self.c1 * self.velocities + 
                           self.c3 * math.sqrt(sys.minv) * R)

    def step(self, config):
        self._momentum_update(sys, 0.5*self.dt)
        self._position_update(sys, 0.5*self.dt)
        self._OU_update(sys, self.dt)
        self._position_update(sys, 0.5*self.dt)
        self._momentum_update(sys, 0.5*self.dt)
