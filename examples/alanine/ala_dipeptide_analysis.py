import sys
import os
sys.path.append(os.path.abspath('../'))
from Simulator import Simulator
from collectivevariable import CV_Function
from snapshot import Snapshot, Configuration
from volume import LambdaVolumePeriodic
from ensemble import EnsembleFactory as ef
from ensemble import (LengthEnsemble, SequentialEnsemble, AllOutXEnsemble,
                      AllInXEnsemble)
from storage import Storage
from trajectory import Trajectory

import mdtraj as md

if __name__ == "__main__":
    storage = Storage(
        filename="trajectory.nc",
        mode='a'
    )

    psi_atoms = [6,8,14,16]
    psi = CV_Function("psi", md.compute_dihedrals, trajdatafmt="mdtraj",
                      indices=[psi_atoms])

    # restore old computed values
    storage.cv.restore(psi)

    for tnum in range(1,storage.trajectory.count()+1):
        traj = storage.trajectory.load(tnum)
        degrees = 180/3.14159 # psi reports in radians; I think in degrees
        psis = [psi(snap) for snap in traj]
        print max(psis)*degrees


