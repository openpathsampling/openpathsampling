from __future__ import absolute_import
import os

from openpathsampling.storage import Storage

import openpathsampling.engines.openmm as peng

try:
    import mdtraj as md
except ImportError:
    md = None

from .test_helpers import data_filename

def setup_package():
    if not md:
        return
    # this should generate the ala_small_traj.nc file which we'll use for
    # everything else
    # NOTE: tests using this must be skipped if there's no mdtraj!
    mdtrajectory = md.load(data_filename("ala_small_traj.pdb"))

    snapshot = peng.snapshot_from_pdb(data_filename("ala_small_traj.pdb"))

    storage = Storage(
        filename=data_filename("ala_small_traj.nc"),
        template=snapshot,
        mode='w'
    )

    mytraj = peng.trajectory_from_mdtraj(mdtrajectory, simple_topology=True)
    storage.trajectories.save(mytraj)

    storage.close()


def teardown_package():
    if os.path.isfile(data_filename("ala_small_traj.nc")):
        os.remove(data_filename("ala_small_traj.nc"))
