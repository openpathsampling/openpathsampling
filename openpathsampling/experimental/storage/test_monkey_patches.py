import openpathsampling as paths
import pytest

from openpathsampling.experimental.storage.monkey_patches import *
from openpathsampling.experimental.storage.collective_variables import \
    CollectiveVariable

@pytest.fixture
def patching_tps_network():
    cv = CollectiveVariable(lambda x: x.xyz[0][0]).named('cv')
    state_A = paths.CVDefinedVolume(cv, 0, 1)
    state_B = paths.CVDefinedVolume(cv, 9, 10)
    tps_network = paths.TPSNetwork(state_A, state_B)
    yield tps_network


def test_double_patch(patching_tps_network):
    orig_dict = patching_tps_network.to_dict()
    global paths
    paths = monkey_patch_all(paths)
    try:
        patched_dict = patching_tps_network.to_dict()
        paths = monkey_patch_all(paths)
        double_patched_dict = patching_tps_network.to_dict()
    finally:
        unpatch(paths)

    assert orig_dict != patched_dict
    assert patched_dict == double_patched_dict


def test_unpatch(patching_tps_network):
    orig_dict = patching_tps_network.to_dict()
    global paths
    paths = monkey_patch_all(paths)
    try:
        patched_dict = patching_tps_network.to_dict()
    finally:
        paths = unpatch(paths)

    unpatched_dict = patching_tps_network.to_dict()
    assert unpatched_dict == orig_dict
    assert unpatched_dict != patched_dict
