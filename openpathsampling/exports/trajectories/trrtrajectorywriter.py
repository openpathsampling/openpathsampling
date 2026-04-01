from openpathsampling.integration_tools import HAS_MDTRAJ
import numpy as np
from .core import TrajectoryWriter

class TRRTrajectoryWriter(TrajectoryWriter):
    """Write a trajectory to a Gromacs TRR file.

    This inludes velocities, so it can be used for restarts.
    """
    def __init__(self):
        if not HAS_MDTRAJ:  # -no-cov-
            raise ImportError("MDTraj is not available")

    @property
    def ext(self):
        return "trr"

    def _write(self, trajectory, filename):
        # this uses some "unofficial" MDTraj API
        import mdtraj as md
        nframes = len(trajectory)
        xyz = np.asarray(trajectory.xyz, dtype=np.float32)
        time = np.asarray([0.0]*nframes, dtype=np.float32)
        box = np.asarray(trajectory.box_vectors, dtype=np.float32)
        lambd = np.asarray([0.0]*nframes, dtype=np.float32)
        vel = np.asarray(trajectory.velocities, dtype=np.float32)
        trr = md.formats.TRRTrajectoryFile(str(filename), 'w')
        step = np.arange(0, nframes, dtype=np.int32)
        trr._write(xyz, time, step, box, lambd, vel)
        trr.close()
