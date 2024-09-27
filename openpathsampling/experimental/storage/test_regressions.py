# miscellaneous regression testing
import pytest
from openpathsampling.tests.test_helpers import data_filename
from openpathsampling.netcdfplus import StorableNamedObject

from openpathsampling.experimental.storage.ops_storage import ops_class_info

class TrajWrapper(StorableNamedObject):
    def __init__(self, traj):
        super(TrajWrapper, self).__init__()
        self.traj = traj


def test_use_mdtraj_codec():
    # smoke test to ensure we *can* use the MDTraj codecs
    md = pytest.importorskip("mdtraj")
    traj = md.load(data_filename("ala_small_traj.pdb"))
    serialized = ops_class_info.default_info.serializer(TrajWrapper(traj))
