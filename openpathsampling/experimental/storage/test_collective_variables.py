import pytest
try:
    from unittest import mock
except ImportError:
    import mock

import numpy as np

import openpathsampling as paths
from openpathsampling.tests.test_helpers import make_1d_traj, data_filename
from openpathsampling.engines import openmm as ops_omm
from openpathsampling.engines.topology import MDTrajTopology

from .collective_variables import *
from ..simstore.test_storable_function import MockBackend
from ..simstore import StorageFunctionHandler
from ..simstore.serialization_helpers import get_uuid


try:
    import mdtraj as md
except ImportError:
    pass
else:
    HAS_MDTRAJ = True


class TestCollectiveVariable(object):
    def setup(self):
        self.cv = CollectiveVariable(lambda x: x.xyz[0][0])
        self.traj = make_1d_traj([1.0, 2.0, 3.0])
        ensemble = paths.LengthEnsemble(3)
        self.sample = paths.Sample(replica=0,
                                   trajectory=self.traj,
                                   ensemble=ensemble)
        self.snap = self.traj[0]

    def test_is_scalar_true(self):
        assert self.cv.is_scalar(self.snap) is True

    @pytest.mark.parametrize('inp_type', ['traj', 'samp'])
    def test_is_scalar_false(self, inp_type):
        inp = {'traj': self.traj, 'samp': self.sample}[inp_type]
        assert self.cv.is_scalar(inp) is False


# TODO: It turns out that trajectories don't currently store their reversed
# partners. This test won't work until that idea is implemented.
# class TestReversibleStorableFunction(object):
    # def setup(self):
        # self.func = ReversibleStorableFunction(
            # func=mock.MagicMock(return_value=3.0),
            # result_type='float'
        # )
        # self.traj = make_1d_traj([1.0, 2.0, 3.0])
        # _ = self.func(self.traj)  # evaluate once; it should be cached
        # # TODO: set up a mock storage, and store results in there
        # backend = MockBackend()
        # self.storage = mock.NonCallableMock(backend=backend)
        # sf_handler = StorageFunctionHandler(self.storage)
        # self.storage._sf_handler = sf_handler
        # sf_handler.register_storable_function(self.func)
        # sf_handler.update_cache(self.func.local_cache)
        # pass

    # def test_get_cached(self):
        # pytest.skip()

    # def test_get_storage(self):
        # # this will first clear the cache on self.func, and then try to load
        # pytest.skip()


class TestCoordinateFunctionCV(object):
    def setup(self):
        # On using side_effect: https://stackoverflow.com/a/16162114
        mock_func = mock.MagicMock(side_effect=lambda s: s.xyz[0][0])
        self.func = CoordinateFunctionCV(mock_func)
        values = [1.0, 2.0, 3.0]
        traj = make_1d_traj(values)
        self.inputs = {'traj': traj.reversed, 'snap': traj[0].reversed}
        self.expected = {'traj': list(reversed(values)), 'snap': values[0]}
        _ = self.func(traj)  # run once to cache initial results
        mock_func.reset_mock()

        backend = MockBackend()
        backend.register_storable_function(get_uuid(self.func), 'float')
        self.storage = mock.NonCallableMock(backend=backend)
        self.storage._sf_handler = StorageFunctionHandler(self.storage)
        self.storage._sf_handler.register_function(self.func)
        self.storage._sf_handler.update_cache(self.func.local_cache)

    @pytest.mark.parametrize('inp_type', ['snap', 'traj'])
    def test_get_cached(self, inp_type):
        inputs = self.inputs[inp_type]
        expected = self.expected[inp_type]
        assert self.func(inputs) == expected
        self.func.func.assert_not_called()
        loop = {'snap': [inputs], 'traj': inputs}[inp_type]
        for item in loop:
            assert self.storage.backend.called_load[get_uuid(item)] == 0

    @pytest.mark.parametrize('inp_type', ['snap', 'traj'])
    def test_get_storage(self, inp_type):
        self.storage.backend.add_storable_function_results(
            table_name=get_uuid(self.func),
            result_dict=self.func.local_cache.result_dict
        )
        self.func.local_cache.clear()
        inputs = self.inputs[inp_type]
        expected = self.expected[inp_type]
        assert self.func(inputs) == expected
        self.func.func.assert_not_called()
        loop = {'snap': [inputs], 'traj': inputs}[inp_type]
        for item in loop:
            assert self.storage.backend.called_load[get_uuid(item)] == 1


class TestMDTrajFunctionCV(object):
    def setup(self):
        if not HAS_MDTRAJ:
            pytest.skip("Unable to import MDTraj")
        pytest.importorskip('simtk.unit')

        self.mdt = md.load(data_filename("ala_small_traj.pdb"))
        top = MDTrajTopology(self.mdt.topology)
        self.traj = ops_omm.tools.trajectory_from_mdtraj(self.mdt)
        self.cv = MDTrajFunctionCV(md.compute_distances, top,
                                   atom_pairs=[[0, 1]])


    @pytest.mark.parametrize('inp_type', ['traj', 'snap'])
    def test_mdtraj_processor(self, inp_type):
        inputs = {'traj': self.traj, 'snap': [self.traj[0]]}[inp_type]
        expected = {'traj': self.mdt.xyz,
                    'snap': np.array([self.mdt.xyz[0]])}[inp_type]
        mdtraj_processor = self.cv.func_config.processor_dict['mdtraj']
        mdtraj_trajectory = mdtraj_processor(inputs)
        np.testing.assert_array_equal(mdtraj_trajectory.xyz, expected)

    def test_eval(self):
        by_traj = self.cv(self.traj)
        by_snap = np.array([self.cv(snap) for snap in self.traj])
        np.testing.assert_array_equal(by_traj, by_snap)
        assert by_traj.shape == (len(self.traj),)
