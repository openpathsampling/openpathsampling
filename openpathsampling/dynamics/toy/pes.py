import numpy as np
from openpathsampling.todict import ops_object

from openpathsampling.util.todict import restores_as_full_object


# The decorator @restores_ allows us to restore the object from a JSON
# string completely and can thus be stored automatically

@ops_object
class Toy_PES(object):
    # For now, we only support additive combinations; maybe someday that can
    # include multiplication, too

    def __add__(self, other):
        return Toy_PES_Add(self, other)

    def __sub__(self, other):
        return Toy_PES_Sub(self, other)

    def kinetic_energy(self, sys):
        v = sys.velocities
        m = sys.mass
        return 0.5*np.dot(m, np.multiply(v,v))

class Toy_PES_Combination(Toy_PES):
    def __init__(self, pes1, pes2, fcn, dfdx_fcn):
        self.pes1 = pes1
        self.pes2 = pes2
        self._fcn = fcn
        self._dfdx_fcn = dfdx_fcn

    def V(self, sys):
        return self._fcn(self.pes1.V(sys), self.pes2.V(sys))

    def dVdx(self, sys):
        return self._dfdx_fcn(self.pes1.dVdx(sys), self.pes2.dVdx(sys))

@ops_object
class Toy_PES_Sub(Toy_PES_Combination):
    def __init__(self, pes1, pes2):
        self.pes1 = pes1
        self.pes2 = pes2
        self._fcn = lambda a, b: a - b
        self._dfdx_fcn = lambda a, b: a - b

@ops_object
class Toy_PES_Add(Toy_PES_Combination):
    def __init__(self, pes1, pes2):
        self.pes1 = pes1
        self.pes2 = pes2
        self._fcn = lambda a, b: a + b
        self._dfdx_fcn = lambda a, b: a + b

@ops_object
class HarmonicOscillator(Toy_PES):
    def __init__(self, A, omega, x0):
        self.A = np.array(A)
        self.omega = np.array(omega)
        self.x0 = np.array(x0)

    def V(self, sys):
        dx = sys.positions - self.x0
        k = self.omega*self.omega*sys.mass
        return 0.5*np.dot(self.A * k, dx * dx)

    def dVdx(self, sys):
        dx = sys.positions - self.x0
        k = self.omega*self.omega*sys.mass
        return self.A*k*dx

@ops_object
class Gaussian(Toy_PES):
    ''' Returns the Gaussian given by A*exp(-\sum_i alpha[i]*(x[i]-x0[i])^2)
    '''
    def __init__(self, A, alpha, x0):
        self.A = A
        self.alpha = np.array(alpha)
        self.x0 = np.array(x0)
        self._local_dVdx = np.zeros(self.x0.size)

    def V(self, sys):
        dx = sys.positions - self.x0
        return self.A*np.exp(-np.dot(self.alpha, np.multiply(dx, dx)))

    def dVdx(self, sys):
        dx = sys.positions - self.x0
        exp_part = self.A*np.exp(-np.dot(self.alpha, np.multiply(dx, dx)))
        for i in range(len(dx)):
            self._local_dVdx[i] = -2*self.alpha[i]*dx[i]*exp_part
        return self._local_dVdx

@ops_object
class OuterWalls(Toy_PES):
    def __init__(self, sigma, x0):
        self.sigma = np.array(sigma)
        self.x0 = np.array(x0)
        self._local_dVdx = np.zeros(self.x0.size)

    def V(self, sys):
        dx = sys.positions - self.x0
        myV = 0.0
        for i in range(len(dx)):
            myV += self.sigma[i]*dx[i]**6
        return myV

    def dVdx(self, sys):
        dx = sys.positions - self.x0
        for i in range(len(dx)):
            self._local_dVdx[i] = 6.0*self.sigma[i]*dx[i]**5
        return self._local_dVdx

@ops_object
class LinearSlope(Toy_PES):
    def __init__(self, m, c):
        self.m = m
        self.c = c
        self._local_dVdx = self.m
        self.dim = len(self.m)

    def V(self, sys):
        return np.dot(self.m, sys.positions) + self.c

    def dVdx(self, sys):
        # this is independent of the position
        return self._local_dVdx
