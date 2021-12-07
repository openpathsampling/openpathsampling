import numpy as np
import openpathsampling as paths
from openpathsampling.tests.test_helpers import make_1d_traj

class TPSSystemFixture(object):
    """Fixture to use with TPS simulation setups.

    This acts as a container for the network and move scheme for any path
    sampling fixture, and also provides the method :meth:`.make_trajectory`
    to make it easy to create trajectory segments.
    """
    def __init__(self, network, scheme):
        self.network = network
        self.scheme = scheme

    @staticmethod
    def make_trajectory(lower, upper):
       """Make a trajectory from `lower` to `upper`.

       This makes a trajectory segment with x-values of frames spread by
       1.  For each trajectory, some random value epsilon is added, such
       that the lowest value is ``lower + epsilon`` and the highest value
       is ``upper + epsilon``.

       Parameters
       ----------
       lower : int
           the lower bound x-value of the trajectory
       upper: int
           the upper bound x-value of the trajectory (NOTE: trajectories will
           cross this value, but not cross ``upper + 0.1``)

       Returns
       -------
       :class:`.Trajectory`
           the trajectory segment
       """
       xvals = np.array(range(lower, upper + 1)) + np.random.random()
       return make_1d_traj(xvals)


class TISSystemFixture(TPSSystemFixture):
    """
    """
    def __init__(self, network, scheme, state_bounds):
        super(TISSystemFixture, self).__init__(network, scheme)
        self.lower, self.upper = state_bounds

    def make_tis_trajectory(self, cv_max, lower_bound=None):
        """Make a TIS trajectory with a given maximum x value.

        Parameters
        ----------
        cv_max : int
            maximum allowed value of x
        lower_bound: int or None
            minumum value of x
        """
        if lower_bound is None:
            lower_bound = self.lower - 1

        increasing = self.make_trajectory(lower_bound, cv_max)

        if cv_max >= self.upper:
            return increasing
        else:
            decreasing = increasing.reversed[1:]
            return increasing + decreasing
