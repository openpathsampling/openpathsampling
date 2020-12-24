import numpy as np

from openpathsampling.netcdfplus import StorableObject


# The decorator @restores_ allows us to restore the object from a JSON
# string completely and can thus be stored automatically

class PES(StorableObject):
    """Abstract base class for toy potential energy surfaces.
    """
    # For now, we only support additive combinations; maybe someday that can
    # include multiplication, too

    def __init__(self):
        super(PES, self).__init__()

    def __add__(self, other):
        return PES_Add(self, other)

    def __sub__(self, other):
        return PES_Sub(self, other)

    def kinetic_energy(self, sys):
        """Default kinetic energy implementation.

        Parameters
        ----------
        sys : :class:`.ToyEngine`
            engine contains its state, including velocities and masses
        """
        v = sys.velocities
        m = sys.mass
        return 0.5*np.dot(m, np.multiply(v, v))


class PES_Combination(PES):
    """Mathematical combination of two potential energy surfaces.

    Abstract base class.

    Parameters
    ----------
    pes1 : :class:`.PES`
        first potential energy surface of the combination
    pes2 : :class:`.PES`
        second potential energy surface of the combination
    fcn : function of two variables
        function to combine the PES energies
    dfdx_fcn : function of two variables
        function to combine the PES (first) derivatives
    """
    def __init__(self, pes1, pes2, fcn, dfdx_fcn):
        super(PES_Combination, self).__init__()
        self.pes1 = pes1
        self.pes2 = pes2
        self._fcn = fcn
        self._dfdx_fcn = dfdx_fcn

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
        return self._fcn(self.pes1.V(sys), self.pes2.V(sys))

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
        return self._dfdx_fcn(self.pes1.dVdx(sys), self.pes2.dVdx(sys))


class PES_Sub(PES_Combination):
    """Difference of two potential energy surfaces; pes1 - pes2

    Parameters
    ----------
    pes1 : :class:`.PES`
        first potential energy surface of the combination
    pes2 : :class:`.PES`
        second potential energy surface of the combination
    """
    def __init__(self, pes1, pes2):
        super(PES_Sub, self).__init__(
            pes1,
            pes2,
            lambda a, b: a - b,
            lambda a, b: a - b
            )


class PES_Add(PES_Combination):
    """Sum of two potential energy surfaces; pes1 + pes 2

    Parameters
    ----------
    pes1 : :class:`.PES`
        first potential energy surface of the combination
    pes2 : :class:`.PES`
        second potential energy surface of the combination
    """
    def __init__(self, pes1, pes2):
        super(PES_Add, self).__init__(
            pes1,
            pes2,
            lambda a, b: a + b,
            lambda a, b: a + b
        )


class HarmonicOscillator(PES):
    r"""Simple harmonic oscillator. Independent in each degree of freedom.

    :math:`V(x) = \sum_i A_i * mass_i * omega_i^2 * (x_i - x0_i)^2`

    Parameters
    ----------
    A : list of float
    omega : list of float
    x0 : list of float
    """
    def __init__(self, A, omega, x0):
        super(HarmonicOscillator, self).__init__()
        self.A = np.array(A)
        self.omega = np.array(omega)
        self.x0 = np.array(x0)

    def __repr__(self):  # pragma: no cover
        repr_str = "HarmonicOscillator({obj.A}, {obj.omega}, {obj.x0})"
        return repr_str.format(obj=self)

    def to_dict(self):
        dct = super(HarmonicOscillator, self).to_dict()
        dct['A'] = dct['A'].tolist()
        dct['omega'] = dct['omega'].tolist()
        dct['x0'] = dct['x0'].tolist()
        return dct

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
        dx = sys.positions - self.x0
        k = self.omega*self.omega*sys.mass
        return 0.5*np.dot(self.A * k, dx * dx)

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
        dx = sys.positions - self.x0
        k = self.omega*self.omega*sys.mass
        return self.A*k*dx


class Gaussian(PES):
    r"""Gaussian given by: :math:`A*exp(-\sum_i alpha[i]*(x[i]-x0[i])^2)`

    Parameters
    ----------
    A : float
        amplitude of the Gaussian
    alpha : list of float
        Gaussian width parameter
    x0 : list of float
        center of the Gaussian
    """
    def __init__(self, A, alpha, x0):
        super(Gaussian, self).__init__()
        self.A = A
        self.alpha = np.array(alpha)
        self.x0 = np.array(x0)
        self._local_dVdx = np.zeros(self.x0.size)

    def to_dict(self):
        dct = super(Gaussian, self).to_dict()
        dct['alpha'] = dct['alpha'].tolist()
        dct['x0'] = dct['x0'].tolist()
        return dct

    def __repr__(self):  # pragma: no cover
        return "Gaussian({o.A}, {o.alpha}, {o.x0})".format(o=self)

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
        dx = sys.positions - self.x0
        return self.A*np.exp(-np.dot(self.alpha, np.multiply(dx, dx)))

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
        dx = sys.positions - self.x0
        exp_part = self.A*np.exp(-np.dot(self.alpha, np.multiply(dx, dx)))
        for i in range(len(dx)):
            self._local_dVdx[i] = -2*self.alpha[i]*dx[i]*exp_part
        return self._local_dVdx


class OuterWalls(PES):
    r"""Creates an x**6 barrier around the system.

    :math:`V(x) = \sum_i sigma_i * (x_i - x0_i)^6`

    Parameters
    ----------
    sigma : list of float
        linear scaling of the potential wall
    x0 : list of float
        center of the potential
    """
    def __init__(self, sigma, x0):
        super(OuterWalls, self).__init__()
        self.sigma = np.array(sigma)
        self.x0 = np.array(x0)
        self._local_dVdx = np.zeros(self.x0.size)

    def to_dict(self):
        dct = super(OuterWalls, self).to_dict()
        dct['x0'] = dct['x0'].tolist()
        dct['sigma'] = dct['sigma'].tolist()
        return dct

    def __repr__(self):  # pragma: no cover
        return "OuterWalls({o.sigma}, {o.x0})".format(o=self)

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
        dx = sys.positions - self.x0
        myV = 0.0
        for i in range(len(dx)):
            myV += self.sigma[i]*dx[i]**6
        return myV

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
        dx = sys.positions - self.x0
        for i in range(len(dx)):
            self._local_dVdx[i] = 6.0*self.sigma[i]*dx[i]**5
        return self._local_dVdx


class LinearSlope(PES):
    r"""Linear potential energy surface.  :math:`V(x) = \sum_i m_i * x_i + c`

    Parameters
    ----------
    m : list of float
        slope
    c : float
        energy offset
    """
    def __init__(self, m, c):
        super(LinearSlope, self).__init__()
        self.m = m
        self.c = c
        self._local_dVdx = self.m
        self.dim = len(self.m)

    def __repr__(self):  # pragma: no cover
        return "LinearSlope({o.m}, {o.c})".format(o=self)

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
        return np.dot(self.m, sys.positions) + self.c

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
        # this is independent of the position
        return self._local_dVdx
