from opentis.storage import Storage, TrajectoryStorage, SnapshotStorage
from opentis.trajectory import Trajectory
from opentis.dynamics_engine import DynamicsEngine
import mdtraj as md

from test_helpers import data_filename

def setup_package():
    # this should generate the trajectory.nc file which we'll use for
    # everything else
    mdtrajectory = md.load(data_filename("ala_small_traj.pdb"))
    engine = DynamicsEngine()
    engine.storage = Storage(
        topology_file=data_filename("ala_small_traj.pdb"),
        filename=data_filename("ala_small_traj.nc"),
        mode='w'
    )
    engine.storage.engine = engine
    Trajectory.storage = engine.storage
    mytraj = Trajectory.from_mdtraj(mdtrajectory)
    engine.storage.trajectory.save(mytraj)


def teardown_package():
    # this should delete the trajectory.nc file
    pass
