"""
Classes to use as fixtures in pytest.

This module includes tools to create containers that are useful as
testsystem-style fixtures. That is, there are connected tools for network,
move scheme, and trajectory generation that may be specific to that network.
"""

import pytest
import numpy as np
import openpathsampling as paths
from openpathsampling.tests.test_helpers import make_1d_traj

DEFAULT_CV = paths.FunctionCV("x", lambda s: s.xyz[0][0])
DEFAULT_ENGINE = paths.engines.NoEngine({'n_spatial': 3, 'n_atoms': 1})


class TPSSystemFixture(object):
    """Fixture to use with TPS simulation setups.

    This acts as a container for the network and move scheme for any path
    sampling fixture, and also provides the method :meth:`.make_trajectory`
    to make it easy to create trajectory segments.

    Parameters
    ----------
    network : :class:`.TransitionNetwork`
        the network associated with this setup
    scheme : :class:`.MoveScheme`
        the move scheme associated with this setup
    """
    def __init__(self, network, scheme):
        self.network = network
        self.scheme = scheme

    @staticmethod
    def make_trajectory(lower, upper):
        """Make a trajectory from `lower` to `upper`.

        This makes a trajectory segment with x-values of frames spread by
        1.  For each trajectory, some random value epsilon is added (with ``0
        <= epsilon < 1``), such that the lowest value is ``lower + epsilon``
        and the highest value is ``upper + epsilon``. This is designed to
        facilitate test debugging, since all frames from the the same
        trajectory will have the same value of epsilon.

        Note that inputs here should be integers. This is very specifically
        for toy models where the boundaries are at integer values (though the
        CVs are allowed to take on arbitrary floats).

        Parameters
        ----------
        lower : int
            the lower bound x-value of the trajectory
        upper: int
            the upper bound x-value of the trajectory (NOTE: trajectories will
            cross this value, but not cross ``upper + 1``)

        Returns
        -------
        :class:`.Trajectory`
            the trajectory segment
        """
        xvals = np.array(range(lower, upper + 1)) + np.random.random()
        return make_1d_traj(xvals)


class TISSystemFixture(TPSSystemFixture):
    """Fixture to use with simple TIS simulation setups.

    This is a container for the network and move scheme for simple two-state
    TIS fixture setups. It provides the :meth:`.make_trajectory` and
    :meth:`.make_tis_trajectory` methods.

    Parameters
    ----------
    network : :class:`.TransitionNetwork`
        the network associated with this setup
    scheme : :class:`.MoveScheme`
        the move scheme associated with this setup
    state_bounds : Tuple[int, int]
        the bounds (lower, upper) that mark the bounds of the transition
        region, i.e., the lower number is the upper bound of the initial
        state, and the upper number is the lower bound of the final state.
    """
    def __init__(self, network, scheme, state_bounds):
        super(TISSystemFixture, self).__init__(network, scheme)
        self.lower, self.upper = state_bounds

    def make_tis_trajectory(self, cv_max, lower_bound=None):
        """Make a TIS trajectory with a given maximum x value.

        As with :method:`.make_trajectory`, this will create a trajectory
        with frames where the x values are separated by 1, with the
        trajectory offset by a random constant epsilon with `0 <= epsilon <
        1`, which allows related frames to by identified when debugging
        analysis tests.

        Parameters
        ----------
        cv_max : int
            maximum allowed value of x; if less than ``self.upper``, then
            this trajectory will return to the lower state, otherwise it
            terminates after reaching ``cv_max``.
        lower_bound: int or None
            minumum value of x; if None, use the ``self.lower - 1``, where
            ``self.lower`` is the lower state bound.
        """
        if lower_bound is None:
            lower_bound = self.lower - 1

        increasing = self.make_trajectory(lower_bound, cv_max)

        if cv_max >= self.upper:
            return increasing
        else:
            decreasing = increasing.reversed[1:]
            return increasing + decreasing


def make_fixture(fixture_type, make_network, scheme_type, **kwargs):
    """Create a pytest fixture for a given setup.

    This will always use ``DEFAULT_ENGINE``, which is a :class:`.NoEngine`
    with 3 spatial dimensions and 1 atom.

    Parameters
    ----------
    fixture_type : type
        the type of fixture (class object) to create; usually either
        :class:`.TPSSystemFixture` or :class:`.TISSystemFixture`
    make_network : Callable[[], :class:`.TransitionNetwork`]
        callable (taking no arguments) that creates that desired transition
        network.
    scheme_type : Callable
        type of move scheme to create, must take arguments ``network`` and
        ``engine`` and return move scheme (usually a scheme type, such as
        :class:`.DefaultScheme` or :class:`.OneWayShootingMoveScheme`).
    kwargs :
        additional keyword arguments are passed to ``fixture_type``
    """
    # TODO: consider whether to only have the pytest.fixture decoration in
    # the conftest.py; that may be better for most use cases
    @pytest.fixture
    def inner():
        network = make_network()
        scheme = scheme_type(network=network, engine=DEFAULT_ENGINE)
        scheme.move_decision_tree()  # to build the movers
        fixture = fixture_type(network=network, scheme=scheme, **kwargs)
        return fixture
    return inner
