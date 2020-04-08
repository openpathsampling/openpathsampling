import pytest
import openpathsampling as paths
from openpathsampling.pathsimulators.adaptive_multilevel_splitting import *
from .test_helpers import make_1d_traj

class TestInterfaceSetParametrizedVolume(object):
    def setup(self):
        paths.InterfaceSet._reset()
        cv = paths.FunctionCV("Id", lambda snap: snap.xyz[0][0])
        self.param_vol = \
                InterfaceSetParametrizedVolume.from_increasing_cv(cv)

    def test_initialization(self):
        assert self.param_vol._volume is None
        assert self.param_vol.cv is not None
        assert self.param_vol.cv_max is not None

    def test_set_parameters(self):
        assert self.param_vol._volume is None
        self.param_vol.set_parameters({'lambda_i': 1.0})
        assert self.param_vol._volume is not None
        traj = make_1d_traj([0.0, 2.0])
        assert self.param_vol(traj[0])
        assert not self.param_vol(traj[1])

    def test_unset_parameters(self):
        assert self.param_vol._volume is None
        self.param_vol.set_parameters({'lambda_i': 1.0})
        assert self.param_vol._volume is not None
        self.param_vol.unset_parameters()
        assert self.param_vol._volume is None

    def test_volume_for_parameters(self):
        assert self.param_vol._volume is None
        vol = self.param_vol.volume_for_parameters({'lambda_i': 1.0})
        traj = make_1d_traj([0.0, 2.0])
        assert vol(traj[0])
        assert not vol(traj[1])

    def test_call_error(self):
        traj = make_1d_traj([0.0, 2.0])
        with pytest.raises(RuntimeError):
            self.param_vol(traj[0])


class TestAMSInitialization(object):
    def setup(self):
        trajectories = [
            make_1d_traj([-0.1, 0.1, -0.11]),
            make_1d_traj([-0.2, 0.2, 0.1, -0.1])
        ]
        self.intialization = AMSInitialization(trajectories)

    def test_serialization_cycle(self):
        dct = self.intialization.to_dict()
        obj = AMSInitialization.from_dict(dct)
        dct2 = obj.to_dict()
        assert dct == dct2
        pass
