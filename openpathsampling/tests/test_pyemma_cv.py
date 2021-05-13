import openpathsampling as paths
import numpy as np
import os
import pytest
from itertools import combinations
from .test_helpers import data_filename
from openpathsampling.engines.openmm.tools import trajectory_from_mdtraj
pyemma = pytest.importorskip("pyemma")
mdtraj = pytest.importorskip("mdtraj")


class TestPyEMMAFeaturizerCV(object):
    def setup(self):
        self.md_trajectory = mdtraj.load(data_filename("ala_small_traj.pdb"))
        self.ops_trajectory = trajectory_from_mdtraj(self.md_trajectory)
        # output of f.pairs(f.select_Backbone()) in the pyemma_generator
        self.pairs = list(combinations([4, 6, 8, 14, 16, 18], 2))

        def pyemma_generator(f):
            f.add_inverse_distances(f.pairs(f.select_Backbone()))

        self.cv = paths.collectivevariable.PyEMMAFeaturizerCV(
            'pyemma',
            pyemma_generator,
            topology=self.ops_trajectory[0].topology
            ).with_diskcache()
        self.storage = paths.Storage('delete.nc', 'w')
        self.fname = self.storage.filename

    def teardown(self):
        fname = self.fname
        if self.storage.isopen():
            self.storage.close()
        if os.path.exists(fname):
            os.unlink(fname)

    def test_proper_result(self):
        res = self.cv(self.ops_trajectory[0])
        # reshape to take out the single frame
        mdres = 1./mdtraj.compute_distances(self.md_trajectory[0],
                                            self.pairs).reshape((-1,))
        assert np.allclose(res, mdres)

    def test_storage_cycle(self):
        # test that pyemma CV can be stored->loaded and still works
        # save a frame
        self.storage.save(self.ops_trajectory[0])
        self.storage.save(self.cv)
        self.storage.close()
        storage = paths.Storage(self.fname, "r")
        # Check if loading by name is still allowed
        cv2 = storage.cvs["pyemma"]
        # Assert that we made a new object
        assert cv2 is not self.cv
        res = self.cv(self.ops_trajectory)
        res2 = cv2(self.ops_trajectory)
        assert np.allclose(res, res2)
        storage.close()

    def test_storage_cycle_results(self):
        # test that pyemma CV results can be stored->loaded and still works
        self.storage.save(self.ops_trajectory[0])
        self.storage.save(self.cv)
        # Generate full results
        _ = self.cv(self.ops_trajectory)
        self.storage.close()
        storage = paths.Storage(self.fname, "r")
        cv2 = storage.cvs["pyemma"]
        # Load in the values
        values = storage.stores["cv%d" % storage.idx(cv2)].variables['value']

        assert values.shape == (0, len(self.pairs))
        assert values.var_type == "numpy.float32"

        feat = storage.vars["attributes_json"][storage.idx(cv2)]
        erg = feat(storage.snapshots)
        res = self.cv(storage.snapshots)
        assert np.allclose(erg, res)
        storage.close()

    def test_scalar(self):
        # Test that PyEmmaCV returns a scalar from a frame
        res = self.cv(self.ops_trajectory[0])
        assert len(res.shape) == 1
        assert res.shape[0] == len(self.pairs)

    def test_iterable(self):
        # Test that PyEmmaCV returns an iteratable from an iteratable
        res = self.cv(self.ops_trajectory)
        mdres = 1./mdtraj.compute_distances(self.md_trajectory, self.pairs)
        assert res.shape == (len(self.ops_trajectory), len(self.pairs))
        assert np.allclose(res, mdres)
