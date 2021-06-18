# extra tests for PathSampling when running SimStore
from openpathsampling.tests import test_pathsimulator
from openpathsampling.experimental.storage import Storage, monkey_patch_all
from openpathsampling.experimental.storage.monkey_patches import unpatch
import openpathsampling as paths

def setup_module():
    import openpathsampling as paths
    paths = monkey_patch_all(paths)

def teardown_module():
    import openpathsampling as paths
    paths = unpatch(paths)

class TestPathSampling(test_pathsimulator.TestPathSampling):
    def test_save_scheme(self, tmpdir):
        filename = tmpdir.join("temp.db")
        storage = Storage(filename, mode='w')
        assert len(storage.schemes) == 0
        sim = paths.PathSampling(storage=storage,
                                 move_scheme=self.scheme,
                                 sample_set=self.init_cond)
        assert len(storage.schemes) == 1

