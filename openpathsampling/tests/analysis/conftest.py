from functools import partial

import pytest
import openpathsampling as paths

from openpathsampling.tests.analysis.fixtures import (
    TISSystemFixture, TPSSystemFixture
)

_DEFAULT_CV = paths.FunctionCV("x", lambda s: s.xyz[0][0])
_DEFAULT_ENGINE = paths.engines.NoEngine({'n_spatial': 3, 'n_atoms': 1})


def unidirectional_tis_network():
    """Fixture for unidirectional TIS with the default (RETIS) scheme.

    This has states defined as :math:`x < 0` and :math:`x \ge 10`. The
    interfaces are at 0, 3, and 6.
    """
    paths.InterfaceSet._reset()
    state_A = paths.CVDefinedVolume(_DEFAULT_CV, float("-inf"), 0)
    state_B = paths.CVDefinedVolume(_DEFAULT_CV, 10, float("inf"))
    interfaces = paths.VolumeInterfaceSet(_DEFAULT_CV, float("-inf"),
                                          [0, 3, 6])
    network = paths.MISTISNetwork([(state_A, interfaces, state_B)])
    return network


def two_state_tps_network():
    state_A = paths.CVDefinedVolume(_DEFAULT_CV, float("-inf"), 0)
    state_B = paths.CVDefinedVolume(_DEFAULT_CV, 10, float("inf"))
    network = paths.TPSNetwork(state_A, state_B)
    return network


def make_scheme(network, scheme_type):
    scheme = scheme_type(network=network, engine=_DEFAULT_ENGINE)
    scheme.move_decision_tree()  # to build the movers
    return scheme


def make_fixture(fixture_type, make_network, scheme_type, **kwargs):
    @pytest.fixture
    def inner():
        network = make_network()
        scheme = make_scheme(network, scheme_type)
        fixture = fixture_type(network=network, scheme=scheme, **kwargs)
        return fixture
    return inner


default_unidirectional_tis = make_fixture(
    fixture_type=TISSystemFixture,
    make_network=unidirectional_tis_network,
    scheme_type=paths.DefaultScheme,
    state_bounds=(0, 10)
)


default_two_state_tps = make_fixture(
    fixture_type=TPSSystemFixture,
    make_network=two_state_tps_network,
    scheme_type=paths.OneWayShootingMoveScheme,
)
