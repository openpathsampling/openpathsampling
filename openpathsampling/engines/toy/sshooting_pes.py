import numpy as np
from .pes import PES

class DoubleWell(PES):
    """Simple double-well potential. Independent in each degree of freedom.
    
    V(x) = \sum_i A_i * (x_i**2 - x0_i**2)**2

    WARNING: Two minima only in one dimension, otherwise there are more!

    Parameters
    ----------
    A : list of float
        potential prefactor for in each degree of freedom.
    x0 : list of float
        minimum position in each degree of freedom.
    """
    def __init__(self, A, x0):
        super(DoubleWell, self).__init__()
        self.A = np.array(A)
        self.x0 = np.array(x0)

    def V(self, sys):
        """Potential energy

        Parameters
        ----------
        sys : :class:`.ToyEngine`
            engine contains its state, including velocities and masses

        Returns
        -------
        float
            the potential energy
        """
        dx2 = sys.positions * sys.positions - self.x0 * self.x0
        return np.dot(self.A, dx2 * dx2)

    def dVdx(self, sys):
        """Derivative of potential energy (-force)

        Parameters
        ----------
        sys : :class:`.ToyEngine`
            engine contains its state, including velocities and masses

        Returns
        -------
        np.array
            the derivatives of the potential at this point
        """
        dx2 = sys.positions * sys.positions - self.x0 * self.x0
        return 4 * self.A * sys.positions * dx2
