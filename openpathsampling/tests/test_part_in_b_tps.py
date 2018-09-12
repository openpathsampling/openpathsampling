from __future__ import absolute_import
from nose.tools import (assert_equal, assert_not_equal, assert_almost_equal,
                        raises)
import logging
from openpathsampling.high_level.part_in_b_tps import *
from .test_network import TestFixedLengthTPSNetwork

logging.getLogger('openpathsampling.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.ensemble').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.netcdfplus').setLevel(logging.CRITICAL)

class TestPartInBNetwork(TestFixedLengthTPSNetwork):
    NetworkType = PartInBFixedLengthTPSNetwork

    def test_allow_self_transitions_false_ABX(self):
        network = self.NetworkType.from_states_all_to_all(
            self.states, allow_self_transitions=False, **self.std_kwargs
        )
        assert_equal(len(network.sampling_ensembles), 1)
        ensemble = network.sampling_ensembles[0]
        assert_equal(ensemble(self.traj['AB0']), True)
        assert_equal(ensemble(self.traj['ABA']), True)

    def test_allow_self_transitions_true_ABX(self):
        network = self.NetworkType.from_states_all_to_all(
            self.states, allow_self_transitions=True, **self.std_kwargs
        )
        assert_equal(len(network.sampling_ensembles), 1)
        ensemble = network.sampling_ensembles[0]
        assert_equal(ensemble(self.traj['AB0']), True)
        assert_equal(ensemble(self.traj['ABA']), True)
