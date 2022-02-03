# TODO: Currently this is located in tests/analysis/conftest.py. However,
# this might be more suitable for a higher-level conftest.py, or perhaps it
# should be moved into `fixtures` subdirectory of tests with the necessary
# objects imported into the main tests/conftest.py, for use across the test
# suite.
import openpathsampling as paths

from openpathsampling.tests.analysis.utils.fixture_classes import (
    TISSystemFixture, TPSSystemFixture, make_fixture, DEFAULT_CV
)


def unidirectional_tis_network():
    r"""Fixture for unidirectional TIS with the default (RETIS) scheme.

    This has states defined as initial state :math:`x < 0` and final state
    :math:`x \ge 10`. The interfaces are at :math:`x=0`, :math:`x=3`, and
    :math:`x=6`.
    """
    paths.InterfaceSet._reset()
    state_A = paths.CVDefinedVolume(DEFAULT_CV, float("-inf"), 0)
    state_B = paths.CVDefinedVolume(DEFAULT_CV, 10, float("inf"))
    interfaces = paths.VolumeInterfaceSet(DEFAULT_CV, float("-inf"),
                                          [0, 3, 6])
    network = paths.MISTISNetwork([(state_A, interfaces, state_B)])
    return network


def two_state_tps_network():
    state_A = paths.CVDefinedVolume(DEFAULT_CV, float("-inf"), 0)
    state_B = paths.CVDefinedVolume(DEFAULT_CV, 10, float("inf"))
    network = paths.TPSNetwork(state_A, state_B)
    return network


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
