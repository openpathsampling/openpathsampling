import os
import numpy as np

from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        assert_almost_equal, raises)
from nose.plugins.skip import Skip, SkipTest
from test_helpers import (
    true_func, assert_equal_array_array, make_1d_traj, data_filename
)

import openpathsampling as paths

import openpathsampling.engines.toy as peng
from openpathsampling.analysis.network import *
from openpathsampling import VolumeFactory as vf

import logging

logging.getLogger('openpathsampling.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.ensemble').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.netcdfplus').setLevel(logging.CRITICAL)


class testMSTISNetwork(object):
    def setup(self):
        xval = paths.CV_Function(name="xA", f=lambda s : s.xyz[0][0])
        self.stateA = paths.CVRangeVolume(xval, float("-inf"), -0.5)
        self.stateB = paths.CVRangeVolume(xval, -0.1, 0.1)
        self.stateC = paths.CVRangeVolume(xval, 0.5, float("inf"))

        ifacesA = vf.CVRangeVolumeSet(xval, float("-inf"), [-0.5, -0.4, -0.3])
        ifacesB = vf.CVRangeVolumeSet(xval, [-0.2, -0.15, -0.1], [0.2, 0.15, 0.1])
        ifacesC = vf.CVRangeVolumeSet(xval, [0.5, 0.4, 0.3], float("inf"))

        self.traj = {}
        self.traj['AA'] = make_1d_traj(
            coordinates=[-0.51, -0.49, -0.52],
            velocities=[1.0]*3
        )
        self.traj['AB'] = make_1d_traj(
            coordinates=[-0.51, -0.25, 0.0],
            velocities=[1.0]*3
        )
        self.traj['BA'] = make_1d_traj(
            coordinates=[0.0, -0.15, -0.35, -0.52],
            velocities=[-1.0]*4
        )
        self.traj['BB'] = make_1d_traj(
            coordinates=[0.0, -0.25, 0.25, 0.02],
            velocities=[1.0]*4
        )
        self.traj['BC'] = make_1d_traj(
            coordinates=[0.01, 0.16, 0.25, 0.53],
            velocities=[1.0]*4
        )
        self.traj['CB'] = make_1d_traj(
            coordinates=[0.52, 0.25, -0.01],
            velocities=[-1.0]*3
        )
        self.traj['CC'] = make_1d_traj(
            coordinates=[0.51, 0.35, 0.55],
            velocities=[1.0]*3
        )
        # A->C magically jumps over B
        self.traj['AC'] = make_1d_traj(
            coordinates=[-0.51, -0.25, 0.25, 0.51],
            velocities=[1.0]*4
        )
        self.traj['CA'] = make_1d_traj(
            coordinates=[0.52, 0.22, -0.22, -0.52],
            velocities=[1.0]*4
        )

        self.mstis = MSTISNetwork([
            (self.stateA, ifacesA, xval),
            (self.stateB, ifacesB, xval),
            (self.stateC, ifacesC, xval)
        ])

    def test_all_states(self):
        assert_equal(set(self.mstis.all_states), 
                     set([self.stateA, self.stateB, self.stateC]))

    def test_trajectories(self):
        # TODO; make this test fully comprehensive? (loop over all
        # possibilities?)
        ensA0 = self.mstis.from_state[self.stateA].ensembles[0]
        ensAm = self.mstis.from_state[self.stateA].ensembles[-1]
        assert_equal(ensA0(self.traj['AA']), True)
        assert_equal(ensAm(self.traj['AA']), False)
        assert_equal(ensAm(self.traj['AB']), True)
        ensB0 = self.mstis.from_state[self.stateB].ensembles[0]
        assert_equal(ensB0(self.traj['AB']), False)
        assert_equal(ensB0(self.traj['BA']), True)
        assert_equal(ensB0(self.traj['BB']), True)
        assert_equal(ensB0(self.traj['BC']), True)
        assert_equal(ensB0(self.traj['AC']), False)
        ensC0 = self.mstis.from_state[self.stateC].ensembles[0]
        assert_equal(ensC0(self.traj['CC']), True)
        assert_equal(ensC0(self.traj['CB']), True)
        assert_equal(ensC0(self.traj['CA']), True)
        assert_equal(ensC0(self.traj['BC']), False)
        assert_equal(ensC0(self.traj['AC']), False)
        assert_equal(ensC0(self.traj['BB']), False)

    def test_ms_outers(self):
        for traj_label in ['AB', 'BA', 'AC', 'CA', 'BC', 'CB']:
            assert_equal(self.mstis.ms_outers[0](self.traj[traj_label]), True)

    def test_sampling_ensembles(self):
        assert_equal(len(self.mstis.from_state[self.stateA].ensembles), 2)
        assert_equal(len(self.mstis.from_state[self.stateB].ensembles), 2)
        assert_equal(len(self.mstis.from_state[self.stateC].ensembles), 2)

    def test_autonaming(self):
        assert_equal(self.stateA.name, "A")
        assert_equal(self.stateB.name, "B")
        assert_equal(self.stateC.name, "C")

        # check that (1) given names stay unchanged; (2) code knows to skip
        # over any default names that have been assigned (i.e., it renames
        # stateC to "C", not to "A"

        # force renaming to weirdness
        self.stateA.name = "B"
        self.stateB.name = "A"
        self.stateC._name = ""
        xval = paths.CV_Function(name="xA", f=lambda s : s.xyz[0][0])
        ifacesA = vf.CVRangeVolumeSet(xval, float("-inf"), [-0.5, -0.4, -0.3])
        ifacesB = vf.CVRangeVolumeSet(xval, [-0.2, -0.15, -0.1], [0.2, 0.15, 0.1])
        ifacesC = vf.CVRangeVolumeSet(xval, [0.5, 0.4, 0.3], float("inf"))
        new_network = MSTISNetwork([
            (self.stateA, ifacesA, xval),
            (self.stateB, ifacesB, xval),
            (self.stateC, ifacesC, xval)
        ])
        assert_equal(self.stateA.name, "B")
        assert_equal(self.stateB.name, "A")
        assert_equal(self.stateC.name, "C")


class testTPSNetwork(object):
    def setup(self):
        from test_helpers import CallIdentity
        xval = CallIdentity()
        self.stateA = paths.CVRangeVolume(xval, float("-inf"), -0.5)
        self.stateB = paths.CVRangeVolume(xval, -0.1, 0.1)
        self.stateC = paths.CVRangeVolume(xval, 0.5, float("inf"))
        self.states = [self.stateA, self.stateB, self.stateC]
        
    # define all the test networks as properties: we can do something
    # similar then for the fixed path length, and just need to override
    # these properties to get all the same tests
    @property
    def network2a(self):
        return TPSNetwork(initial_states=[self.stateA], 
                          final_states=[self.stateB])

    @property
    def network2b(self):
        return TPSNetwork(initial_states=self.stateA,
                          final_states=self.stateB)

    @property
    def network2c(self):
        return TPSNetwork.from_state_pairs([(self.stateA, self.stateB)])

    @property
    def network3a(self):
        return TPSNetwork(initial_states=self.states,
                          final_states=self.states)

    @property
    def network3b(self):
        return TPSNetwork.from_states_all_to_all(self.states)

    @property
    def network3c(self):
        return TPSNetwork.from_state_pairs([
            (self.stateA, self.stateB), (self.stateA, self.stateC),
            (self.stateB, self.stateA), (self.stateB, self.stateC),
            (self.stateC, self.stateA), (self.stateC, self.stateB)
        ])

    def test_initialization_2state(self):
        network2a = self.network2a
        assert_equal(len(network2a.sampling_transitions), 1)
        assert_equal(len(network2a.transitions), 1)
        network2b = self.network2b
        assert_equal(len(network2b.sampling_transitions), 1)
        assert_equal(len(network2b.transitions), 1)
        network2c = self.network2c
        assert_equal(len(network2c.sampling_transitions), 1)
        assert_equal(len(network2c.transitions), 1)
        assert_equal(set(network2a.all_states), set(network2b.all_states))
        assert_equal(set(network2b.all_states), set(network2c.all_states))
        assert_equal(set(network2a.all_states), 
                     set([self.stateA, self.stateB]))


    def test_initialization_3state(self):
        network3a = self.network3a
        assert_equal(len(network3a.sampling_transitions), 1)
        assert_equal(len(network3a.transitions), 6)
        network3b = self.network3b
        assert_equal(len(network3b.sampling_transitions), 1)
        assert_equal(len(network3b.transitions), 6)
        network3c = self.network3c
        assert_equal(len(network3c.sampling_transitions), 1)
        assert_equal(len(network3c.transitions), 6)

    def test_storage(self):
        import os
        fname = data_filename("tps_network_storage_test.nc")
        if os.path.isfile(fname):
            os.remove(fname)
        topol = peng.Topology(n_spatial=1, masses=[1.0], pes=None)
        self.template = peng.Snapshot(coordinates=np.array([[0.0]]),
                                       velocities=np.array([[0.0]]),
                                       topology=topol)
        states = [self.stateA, self.stateB, self.stateC]
        network_a = TPSNetwork(initial_states=states, final_states=states)
        assert_equal(len(network_a.sampling_transitions), 1)
        assert_equal(len(network_a.transitions), 6)
        storage_w = paths.storage.Storage(fname, "w", self.template)
        storage_w.save(network_a)
        storage_w.sync_all()

        storage_r = paths.storage.AnalysisStorage(fname)
        network_b = storage_r.networks[0]
        assert_equal(len(network_b.sampling_transitions), 1)
        assert_equal(len(network_b.transitions), 6)

        if os.path.isfile(fname):
            os.remove(fname)

class testFixedLengthTPSNetwork(testTPSNetwork):
    @property
    def network2a(self):
        return FixedLengthTPSNetwork(initial_states=[self.stateA],
                                     final_states=[self.stateB],
                                     length=10)

    @property
    def network2b(self):
        return FixedLengthTPSNetwork(initial_states=self.stateA,
                                     final_states=self.stateB,
                                     length=10)

    @property
    def network2c(self):
        return FixedLengthTPSNetwork.from_state_pairs(
            [(self.stateA, self.stateB)], 
            length=10
        )

    @property
    def network3a(self):
        return FixedLengthTPSNetwork(initial_states=self.states,
                                     final_states=self.states,
                                     length=10)

    @property
    def network3b(self):
        return FixedLengthTPSNetwork.from_states_all_to_all(self.states,
                                                            length=10)

    @property
    def network3c(self):
        return FixedLengthTPSNetwork.from_state_pairs(
            [(self.stateA, self.stateB), (self.stateA, self.stateC),
             (self.stateB, self.stateA), (self.stateB, self.stateC),
             (self.stateC, self.stateA), (self.stateC, self.stateB)],
            length=10
        )

    def test_lengths(self):
        for network in [self.network2a, self.network2b, self.network2c,
                        self.network3a, self.network3b, self.network3c]:
            assert_equal(network.sampling_transitions[0].length, 10)
            assert_equal(network.transitions.values()[0].length, 10)
