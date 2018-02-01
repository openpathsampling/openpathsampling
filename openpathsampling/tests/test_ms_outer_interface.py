from __future__ import absolute_import
from builtins import object
from nose.tools import (assert_equal, assert_not_equal, assert_almost_equal,
                        raises)
from nose.plugins.skip import Skip, SkipTest
from .test_helpers import (
    true_func, assert_equal_array_array, make_1d_traj, data_filename,
    assert_items_equal
)

import openpathsampling as paths
from openpathsampling.high_level.interface_set import GenericVolumeInterfaceSet

import logging
logging.getLogger('openpathsampling.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.ensemble').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.netcdfplus').setLevel(logging.CRITICAL)

class TestMSOuterTISInterface(object):
    def setup(self):
        paths.InterfaceSet._reset()
        self.cv_inc = paths.FunctionCV(name="inc", f=lambda s: s.xyz[0][0])
        self.cv_dec = paths.FunctionCV(name="dec", 
                                        f=lambda s: 1.0-s.xyz[0][0])
        self.lambdas = [0.0, 0.1, 0.2, 0.3]
        self.interfaces_inc = paths.VolumeInterfaceSet(cv=self.cv_inc,
                                                       minvals=float("-inf"),
                                                       maxvals=self.lambdas)
        self.interfaces_dec = paths.VolumeInterfaceSet(cv=self.cv_dec,
                                                       minvals=float("-inf"),
                                                       maxvals=self.lambdas)
        self.stateA = paths.CVDefinedVolume(self.cv_inc, float("-inf"), 0.0
                                         ).named("A")
        self.stateB = paths.CVDefinedVolume(self.cv_dec, float("-inf"), 0.0
                                         ).named("B")
        self.network = paths.MISTISNetwork([
            (self.stateA, self.interfaces_inc, self.stateB),
            (self.stateB, self.interfaces_dec, self.stateA)
        ])
        self.volumes = [self.interfaces_inc.new_interface(0.5),
                        self.interfaces_dec.new_interface(0.4)]

        self.ms_outer_explicit = paths.MSOuterTISInterface(
            interface_sets=[self.interfaces_inc, self.interfaces_dec],
            volumes=self.volumes,
            lambdas=[0.5, 0.4]
        )

        self.ms_outer = paths.MSOuterTISInterface.from_lambdas(
            {self.interfaces_inc: 0.5, self.interfaces_dec: 0.4}
        )

        # TODO: temporary hack until networks working; remove after
        self.post_network = paths.MSOuterTISInterface.from_lambdas(
            {t.interfaces: {self.cv_inc: 0.5,
                            self.cv_dec: 0.4}[t.interfaces.cv]
             for t in self.network.sampling_transitions}
        )

    def test_initialization(self):
        by_lambda = self.ms_outer
        explicit = self.ms_outer_explicit
        for iface_set in explicit.interface_sets:
            assert_equal(by_lambda.volume_for_interface_set(iface_set),
                         explicit.volume_for_interface_set(iface_set))
            assert_equal(by_lambda.lambda_for_interface_set(iface_set),
                         explicit.lambda_for_interface_set(iface_set))

        assert_equal(len(explicit.volumes), len(by_lambda.volumes))
        assert_equal(len(explicit.lambdas), len(by_lambda.lambdas))
        assert_items_equal(set(explicit.interface_sets),
                           set(by_lambda.interface_sets))


    def test_volume_for_interface_set(self):
        assert_equal(
            self.ms_outer.volume_for_interface_set(self.interfaces_inc),
            self.volumes[0]
        )
        assert_equal(
            self.ms_outer.volume_for_interface_set(self.interfaces_dec),
            self.volumes[1]
        )

    def test_lambda_for_interface_set(self):
        assert_equal(
            self.ms_outer.lambda_for_interface_set(self.interfaces_inc),
            0.5
        )
        assert_equal(
            self.ms_outer.lambda_for_interface_set(self.interfaces_dec),
            0.4
        )

    def test_relevant_transitions(self):
        extra_set = paths.VolumeInterfaceSet(self.cv_inc, 0.0, [0.2, 0.3])
        extra = paths.TISTransition(self.stateA, self.stateB, extra_set,
                                    self.cv_inc, "fake")
        transitions = self.network.sampling_transitions + [extra]
        # TODO: switch once network is working
        relevant = self.post_network.relevant_transitions(transitions)
        #relevant = self.ms_outer.relevant_transitions(transitions)
        assert_equal(len(relevant), 2)
        assert_equal(set(self.network.sampling_transitions), set(relevant))

    def test_make_ensemble(self):
        transitions = self.network.sampling_transitions
        # TODO: switch once network is working
        ensemble = self.post_network.make_ensemble(transitions)
        #ensemble = self.ms_outer.make_ensemble(transitions)
        test_AA = make_1d_traj([-0.1, 0.2, -0.2])
        test_AXA = make_1d_traj([-0.1, 0.6, -0.2])
        test_BB = make_1d_traj([1.1, 0.9, 1.2])
        test_BXB = make_1d_traj([1.1, 0.5, 1.2])
        test_AXB = make_1d_traj([-0.1, 0.6, 1.1])
        test_BXA = make_1d_traj([1.1, 0.5, -0.1])
        assert_equal(ensemble(test_AA), False)
        assert_equal(ensemble(test_AXA), True)
        assert_equal(ensemble(test_BB), False)
        assert_equal(ensemble(test_BXB), True)
        assert_equal(ensemble(test_BXA), True)
        assert_equal(ensemble(test_AXB), True)

    def test_make_ensemble_with_forbidden(self):
        forbidden = paths.CVDefinedVolume(self.cv_inc, 0.55, 0.65)
        transitions = self.network.sampling_transitions
        # TODO: switch once network is working
        ensemble = self.post_network.make_ensemble(transitions, forbidden)
        #ensemble = self.ms_outer.make_ensemble(transitions, forbidden)
        test_AA = make_1d_traj([-0.1, 0.2, -0.2])
        test_AXA = make_1d_traj([-0.1, 0.7, -0.2])
        test_AFA = make_1d_traj([-0.1, 0.6, -0.2])
        test_BB = make_1d_traj([1.1, 0.9, 1.2])
        test_BXB = make_1d_traj([1.1, 0.5, 1.2])
        test_BFB = make_1d_traj([1.1, 0.6, 1.2])
        test_AXB = make_1d_traj([-0.1, 0.7, 1.1])
        test_AFB = make_1d_traj([-0.1, 0.6, 1.1])
        test_BXA = make_1d_traj([1.1, 0.5, -0.1])
        test_BFA = make_1d_traj([1.1, 0.6, -0.1])
        assert_equal(ensemble(test_AA), False)
        assert_equal(ensemble(test_AXA), True)
        assert_equal(ensemble(test_BB), False)
        assert_equal(ensemble(test_BXB), True)
        assert_equal(ensemble(test_BXA), True)
        assert_equal(ensemble(test_AXB), True)
        assert_equal(ensemble(test_AFA), False)
        assert_equal(ensemble(test_BFB), False)
        assert_equal(ensemble(test_AFB), False)
        assert_equal(ensemble(test_BFA), False)
