from __future__ import absolute_import
from builtins import zip
from builtins import range
from builtins import object

import pytest
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

class TestInterfaceSet(object):
    def setup_method(self):
        paths.InterfaceSet._reset()
        self.cv = paths.FunctionCV(name="x", f=lambda s: s.xyz[0][0])
        self.lambdas = [0.0, 0.1, 0.2, 0.3]
        min_vals= [float("-inf")] * len(self.lambdas)
        self.volumes = [paths.CVDefinedVolume(self.cv, min_v, max_v)
                        for min_v, max_v in zip(min_vals, self.lambdas)]
        self.interface_set = paths.InterfaceSet(self.volumes, self.cv,
                                                self.lambdas)
        self.decreasing = paths.InterfaceSet(list(reversed(self.volumes)),
                                             self.cv,
                                             list(reversed(self.lambdas)))
        self.no_lambda_set = paths.InterfaceSet(self.volumes, self.cv)

    def test_direction(self):
        assert self.interface_set.direction == 1
        assert self.no_lambda_set.direction == 0
        assert self.decreasing.direction == -1

    def test_get_lambda(self):
        for (v, l) in zip(self.volumes, self.lambdas):
            assert self.interface_set.get_lambda(v) == l
            assert self.no_lambda_set.get_lambda(v) is None

    def test_list_behavior(self):
        # len
        assert len(self.interface_set) == 4
        assert len(self.no_lambda_set) == 4
        # getitem, contains
        for i in range(4):
            assert self.volumes[i] == self.interface_set[i]
            assert self.volumes[i] in self.interface_set
        # getitem for slices
        sliced = self.interface_set[0:2]
        for vol in sliced:
            assert (sliced.get_lambda(vol)
                    == self.interface_set.get_lambda(vol))
        # special case of -1 needs to work (used frequently!)
        assert self.volumes[-1] == self.interface_set[-1]
        # iter
        for vol in self.interface_set:
            assert vol in self.volumes
        # reversed
        i = 0
        for vol in reversed(self.interface_set):
            assert vol == self.volumes[3-i]
            i += 1

    def test_no_direction_possible(self):
        min_vals=[-0.1, -0.2, -0.3]
        max_vals=[0.1, 0.2, 0.3]
        volumes = [paths.CVDefinedVolume(self.cv, min_v, max_v)
                   for min_v, max_v in zip(min_vals, max_vals)]
        ifaces = paths.InterfaceSet(volumes)
        assert ifaces.cv is None
        assert ifaces.cv_max is None
        assert ifaces.direction == 0


class TestGenericVolumeInterfaceSet(object):
    def test_sanitize_input(self):
        # this is just to make the rest a little more readable
        sanitize = GenericVolumeInterfaceSet._sanitize_input
        assert (([float("-inf")]*3, [0.0, 0.1, 0.2], 1) \
                 == sanitize(float("-inf"), [0.0, 0.1, 0.2]))
        assert (([0.2, 0.1, 0.0], [float("inf")]*3, -1) \
                 == sanitize([0.2, 0.1, 0.0], float("inf")))
        assert (([-0.1, -0.2], [0.1, 0.2], 0) \
                 == sanitize([-0.1, -0.2], [0.1, 0.2]))
        assert (([0.0, 0.0], [0.1, 0.2], 1) \
                 == sanitize([0.0, 0.0], [0.1, 0.2]))
        assert (([-0.1, -0.2], [0.0, 0.0], -1) \
                 == sanitize([-0.1, -0.2], [0.0, 0.0]))
        # and the idiot case:
        assert (([-0.1, -0.1], [0.1, 0.1], 0) \
                 == sanitize([-0.1, -0.1], [0.1, 0.1]))

    def test_bad_sanitize(self):
        with pytest.raises(RuntimeError):
            GenericVolumeInterfaceSet._sanitize_input([0.0, -0.1],
                                                      [0.1, 0.2, 0.3])


class TestVolumeInterfaceSet(object):
    def setup_method(self):
        paths.InterfaceSet._reset()
        self.cv = paths.FunctionCV(name="x", f=lambda s: s.xyz[0][0])
        self.increasing_set = paths.VolumeInterfaceSet(cv=self.cv,
                                                       minvals=float("-inf"),
                                                       maxvals=[0.0, 0.1])
        self.decreasing_set = paths.VolumeInterfaceSet(cv=self.cv,
                                                       minvals=[0.0, -0.1],
                                                       maxvals=float("inf"))
        self.weird_set = paths.VolumeInterfaceSet(cv=self.cv,
                                                  minvals=[-0.1, -0.2],
                                                  maxvals=[0.1, 0.2])

    def test_initialization(self):
        assert len(paths.InterfaceSet._cv_max_dict) == 1
        cv_max = list(paths.InterfaceSet._cv_max_dict.values())[0]

        assert len(self.increasing_set) == 2
        assert self.increasing_set.direction == 1
        assert self.increasing_set.lambdas == [0.0, 0.1]
        assert self.increasing_set.cv_max == cv_max

        assert len(self.decreasing_set) == 2
        assert self.decreasing_set.direction == -1
        assert self.decreasing_set.lambdas == [0.0, -0.1]
        # TODO: decide what to do about cv_max for decreasing/weird

        assert len(self.weird_set) == 2
        assert self.weird_set.direction == 0
        assert self.weird_set.lambdas is None

    def test_new_interface(self):
        new_iface = self.increasing_set.new_interface(0.25)
        expected = paths.CVDefinedVolume(self.cv, float("-inf"), 0.25)
        assert expected == new_iface

    def test_bad_new_interface(self):
        with pytest.raises(TypeError):
            self.weird_set.new_interface(0.25)

    def test_storage(self):
        import os
        fname = data_filename("interface_set_storage_test.nc")
        if os.path.isfile(fname):
            os.remove(fname)
        template_traj = make_1d_traj([0.0])
        storage_w = paths.Storage(fname, "w")
        storage_w.save(template_traj)
        storage_w.save(self.increasing_set)
        storage_w.sync_all()
        storage_w.close()

        storage_r = paths.AnalysisStorage(fname)
        reloaded = storage_r.interfacesets[0]

        assert_items_equal(reloaded.lambdas, self.increasing_set.lambdas)
        for (truth, beauty) in zip(self.increasing_set, reloaded):
            assert truth == beauty

        for (v, l) in zip(reloaded.volumes, reloaded.lambdas):
            assert reloaded.get_lambda(v) == l

        if os.path.isfile(fname):
            os.remove(fname)


class TestPeriodicVolumeInterfaceSet(object):
    def setup_method(self):
        paths.InterfaceSet._reset()
        self.cv = paths.FunctionCV(name="x", f=lambda s: s.xyz[0][0])
        self.increasing_set = paths.PeriodicVolumeInterfaceSet(
            cv=self.cv,
            minvals=0.0,
            maxvals=[100, 150, 200-360],
            period_min=-180,
            period_max=180
        )

    def test_initialization(self):
        assert self.increasing_set.direction == 1
        assert len(self.increasing_set) == 3
        assert self.increasing_set.lambdas == [100, 150, -160]

    def test_new_interface(self):
        new_iface = self.increasing_set.new_interface(-140)
        expected = paths.PeriodicCVDefinedVolume(self.cv, 0.0, -140, -180, 180)
        assert new_iface == expected

    def test_storage(self):
        import os
        fname = data_filename("interface_set_storage_test.nc")
        if os.path.isfile(fname):
            os.remove(fname)
        template_traj = make_1d_traj([0.0])
        template = template_traj[0]
        storage_w = paths.Storage(fname, "w")
        storage_w.save(template_traj)
        storage_w.save(self.increasing_set)
        storage_w.sync_all()

        storage_r = paths.AnalysisStorage(fname)
        reloaded = storage_r.interfacesets[0]

        assert_items_equal(reloaded.lambdas, self.increasing_set.lambdas)
        assert reloaded.period_min == self.increasing_set.period_min
        assert reloaded.period_max == self.increasing_set.period_max
        for (truth, beauty) in zip(self.increasing_set, reloaded):
            assert truth == beauty

        for (v, l) in zip(reloaded.volumes, reloaded.lambdas):
            assert reloaded.get_lambda(v) == l

        storage_r.close()
        storage_w.close()

        if os.path.isfile(fname):
            os.remove(fname)
