import numpy as np

from openpathsampling.netcdfplus import StorableObject


# The decorator @restores_ allows us to restore the object from a JSON
# string completely and can thus be stored automatically

class PES(StorableObject):
    # For now, we only support additive combinations; maybe someday that can
    # include multiplication, too

    def __init__(self):
        super(PES, self).__init__()

    def __add__(self, other):
        return PES_Add(self, other)

    def __sub__(self, other):
        return PES_Sub(self, other)

    def kinetic_energy(self, sys):
        v = sys.velocities
        m = sys.mass
        return 0.5*np.dot(m, np.multiply(v,v))

class PES_Combination(PES):
    def __init__(self, pes1, pes2, fcn, dfdx_fcn):
        super(PES_Combination, self).__init__()
        self.pes1 = pes1
        self.pes2 = pes2
        self._fcn = fcn
        self._dfdx_fcn = dfdx_fcn

    def V(self, sys):
        return self._fcn(self.pes1.V(sys), self.pes2.V(sys))

    def dVdx(self, sys):
        return self._dfdx_fcn(self.pes1.dVdx(sys), self.pes2.dVdx(sys))

class PES_Sub(PES_Combination):
    def __init__(self, pes1, pes2):
        super(PES_Sub, self).__init__(
            pes1,
            pes2,
            lambda a, b: a - b,
            lambda a, b: a - b
            )

class PES_Add(PES_Combination):
    def __init__(self, pes1, pes2):
        super(PES_Add, self).__init__(
            pes1,
            pes2,
            lambda a, b: a + b,
            lambda a, b: a + b
        )

class HarmonicOscillator(PES):
    def __init__(self, A, omega, x0):
        super(HarmonicOscillator, self).__init__()
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

class Gaussian(PES):
    ''' Returns the Gaussian given by A*exp(-\sum_i alpha[i]*(x[i]-x0[i])^2)
    '''
    def __init__(self, A, alpha, x0):
        super(Gaussian, self).__init__()
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

class OuterWalls(PES):
    def __init__(self, sigma, x0):
        super(OuterWalls, self).__init__()
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

class LinearSlope(PES):
    def __init__(self, m, c):
        super(LinearSlope, self).__init__()
        self.m = m
        self.c = c
        self._local_dVdx = self.m
        self.dim = len(self.m)

    def V(self, sys):
        return np.dot(self.m, sys.positions) + self.c

    def dVdx(self, sys):
        # this is independent of the position
        return self._local_dVdx
