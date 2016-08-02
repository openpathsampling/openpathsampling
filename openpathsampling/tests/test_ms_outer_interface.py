from nose.tools import (assert_equal, assert_not_equal, assert_items_equal,
                        assert_almost_equal, raises)
from nose.plugins.skip import Skip, SkipTest
from test_helpers import (
    true_func, assert_equal_array_array, make_1d_traj, data_filename
)

import openpathsampling as paths
from openpathsampling.high_level.interface_set import GenericVolumeInterfaceSet

import logging
logging.getLogger('openpathsampling.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.ensemble').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.netcdfplus').setLevel(logging.CRITICAL)

class testMSOuterTISInterface(object):
    def setup(self):
        self.cv_inc = paths.CV_Function(name="inc", f=lambda s: s.xyz[0][0])
        self.cv_dec = paths.CV_Function(name="dec", 
                                        f=lambda s: 1.0-s.xyz[0][0])
        self.lambdas = [0.0, 0.1, 0.2, 0.3]
        self.interfaces_inc = paths.VolumeInterfaceSet(cv=self.cv_inc,
                                                       minvals=float("-inf"),
                                                       maxvals=self.lambdas)
        self.interfaces_dec = paths.VolumeInterfaceSet(cv=self.cv_dec,
                                                       minvals=float("inf"),
                                                       maxvals=self.lambdas)
        self.stateA = paths.CVRangeVolume(self.cv_inc, float("-inf"), 0.0
                                         ).named("A")
        self.stateB = paths.CVRangeVolume(self.cv_dec, float("-inf"), 0.0
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

    def test_initialization(self):
        by_lambda = self.ms_outer
        explicit = self.ms_outer_explicit
        for iface_set in explicit.interface_sets:
            assert_equal(by_lambda.volume_for_interface_set(iface_set),
                         explicit.volume_for_interface_set(iface_set))
            assert_equal(by_lambda.lambda_for_interface_set(iface_set),
                         explicit.lambda_for_interface_set(iface_set))

        for (v1, v2) in zip(explicit.volumes, by_lambda.volumes):
            assert_equal(v1, v2)

        for (l1, l2) in zip(explicit.lambdas, by_lambda.lambdas):
            assert_equal(l1, l2)

        assert_items_equal(explicit.interface_sets,
                           by_lambda.interface_sets)


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
        # add an extra transition
        raise SkipTest

    def test_make_ensemble(self):
        raise SkipTest
