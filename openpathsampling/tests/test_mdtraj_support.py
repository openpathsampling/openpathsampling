import logging

from nose.tools import (
    assert_equal, assert_not_equal, assert_items_equal, raises
)
from nose.plugins.skip import SkipTest
import numpy.testing as nptest
from test_helpers import data_filename

import openpathsampling as paths
import mdtraj as md

from openpathsampling.engines.openmm.tools import (
    trajectory_from_mdtraj, trajectory_to_mdtraj
)

logging.getLogger('opentis.trajectory').setLevel(logging.DEBUG)
logging.getLogger('opentis.initialization').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.storage').setLevel(logging.CRITICAL)
logging.getLogger('openpathsampling.netcdfplus').setLevel(logging.CRITICAL)

# Includes tests of trajectory conversion to and from MDTraj; should also
# add tests for topology conversions as well

class testMDTrajSupport(object):
    def setUp(self):
        self.md_trajectory = md.load(data_filename("ala_small_traj.pdb"))

    def test_trajectory_mdtraj_round_trip(self):
        # first test that the ops trajectory works
        ops_trajectory = trajectory_from_mdtraj(self.md_trajectory)
        nptest.assert_allclose(ops_trajectory.xyz, self.md_trajectory.xyz)
        nptest.assert_allclose(ops_trajectory.box_vectors,
                               self.md_trajectory.unitcell_vectors)
    
        # switch back to mdtraj
        md_trajectory_2 = trajectory_to_mdtraj(ops_trajectory,
                                               ops_trajectory.topology.mdtraj)
        nptest.assert_allclose(self.md_trajectory.xyz, md_trajectory_2.xyz)
        nptest.assert_allclose(self.md_trajectory.unitcell_vectors,
                               md_trajectory_2.unitcell_vectors)
        




