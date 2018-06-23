import numpy as np
from .integrators import ToyIntegrator

class OverdampedLangevinIntegrator(ToyIntegrator):
    """Over-damped Langevin integrator

    Simple integrator for over-damped Langevin dynamics.

    Parameters
    ----------
    dt : float
        time step
    temperature : float
        temperature
    D : float
        diffusion constant
    """
    
    dd = None

    def __init__(self, dt, temperature, D):
        super(OverdampedLangevinIntegrator, self).__init__()
        self.dt = dt
        self.temperature = temperature
        self.beta = 1.0 / temperature
        self.D = D
        self.A = D * dt / temperature
        self.R = np.sqrt(2.0 * D * dt)

    def _position_update(self, sys, mydt):
        sys.positions += - self.A * np.array(sys.pes.dVdx(sys)) \
                         + self.R * np.random.normal(size=len(sys.positions))

    def step(self, sys):
        """
        Take an MD step. Update in-place.

        Parameters
        ----------
        sys : :class:`.ToyEngine`
            engine contains its state, including velocities and masses
        """
        self._position_update(sys, self.dt)

