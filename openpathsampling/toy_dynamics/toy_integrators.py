import math

import numpy as np

from openpathsampling.base import StorableNamedObject


class ToyIntegrator(StorableNamedObject):
    def __init__(self):
        super(StorableNamedObject, self).__init__()

class LeapfrogVerletIntegrator(ToyIntegrator):
    """Leapfrog Integrator. Not for actual use; but the momentum and
    position update functions are used in other integrators, so we inherit
    from this.
    """
    
    dd = None

    def __init__(self, dt):
        super(LeapfrogVerletIntegrator, self).__init__()
        self.dt = dt

    def _momentum_update(self, sys, mydt):
        sys.velocities -= sys.pes.dVdx(sys)*sys._minv*mydt

    def _position_update(self, sys, mydt):
        sys.positions += sys.velocities * mydt

    def step(self, sys, nsteps):
        self._position_update(sys, 0.5*self.dt)
        self._momentum_update(sys, self.dt)
        self._position_update(sys, 0.5*self.dt)

class LangevinBAOABIntegrator(LeapfrogVerletIntegrator):
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
        super(LangevinBAOABIntegrator, self).__init__(dt)
        self._beta = 1.0 / temperature
        self._c1 = math.exp(-gamma*dt)
        self._c2 = (1.0-self._c1)/gamma
        self._c3 = math.sqrt((1.0 - self._c1*self._c1) / self._beta)

        self.temperature = temperature
        self.gamma = gamma

    @property
    def beta(self):
        return self._beta

    @property
    def c1(self):
        return self._c1

    @property
    def c2(self):
        return self._c2

    @property
    def c3(self):
        return self._c3


    def _OU_update(self, sys, mydt):
        R = np.random.normal(size=len(sys.velocities))
        sys.velocities = (self._c1 * sys.velocities +
                          self._c3 * np.sqrt(sys._minv) * R)

    def step(self, sys, nsteps):
        self._momentum_update(sys, 0.5*self.dt)
        self._position_update(sys, 0.5*self.dt)
        self._OU_update(sys, self.dt)
        self._position_update(sys, 0.5*self.dt)
        self._momentum_update(sys, 0.5*self.dt)
