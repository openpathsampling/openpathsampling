import numpy as np
from math import exp


class Toy_PES(object):
    # For now, we only support additive combinations; maybe someday that can
    # include multiplication, too
    def __add__(self, other):
        pass

    def __sub__(self, other):
        pass

    def kinetic_energy(self, sys):
        v = sys.velocities
        m = sys.mass
        return np.dot(m, np.multiply(v,v))

class Toy_PES_Combination(Toy_PES):
    def __init__(self, pes1, pes2, fcn):
        pass

    def V(self, sys):
        pass

    def dVdx(self, sys):
        pass


class Gaussian(Toy_PES):
    ''' Returns the Gaussian given by A*exp(-\sum_i alpha[i]*(x[i]-x0[i])^2)
    '''
    def __init__(self, A, alpha, x0):
        self.A = A
        self.alpha = np.array(alpha)
        self.x0 = np.array(x0)
        self.local_dVdx = np.zeros(self.x0.size)

    def V(self, sys):
        dx = sys.positions - self.x0
        return self.A*np.exp(-np.dot(self.alpha, np.multiply(dx, dx)))

    def dVdx(self, sys):
        dx = sys.positions - self.x0
        exp_part = self.A*np.exp(-np.dot(self.alpha, np.multiply(dx, dx)))
        for i in range(len(dx)):
            self.local_dVdx[i] = -2*self.alpha[i]*dx[i]*exp_part
        return self.local_dVdx

class OuterWalls(Toy_PES):
    def __init__(self, sigma, x0):
        self.sigma = np.array(sigma)
        self.x0 = np.array(x0)
        self.local_dVdx = np.zeros(self.x0.size)

    def V(self, sys):
        dx = sys.positions - self.x0
        myV = 0.0
        for i in range(len(dx)):
            myV += self.sigma[i]*dx[i]**6
        return myV

    def dVdx(self, sys):
        dx = sys.positions - self.x0
        for i in range(len(dx)):
            self.local_dVdx[i] = 6.0*self.sigma[i]*dx[i]**5
        return self.local_dVdx


class LinearSlope(Toy_PES):
    def V(self, sys):
        pass

    def dVdx(self, sys):
        pass
