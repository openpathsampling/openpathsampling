def missing_openmm(*args, **kwargs):
    raise RuntimeError("Install OpenMM to use this feature")

try:
    import simtk.openmm
    import simtk.openmm.app
except ImportError:
    HAS_OPENMM = False
    Engine = missing_openmm
    empty_snapshot_from_openmm_topology = missing_openmm
    snapshot_from_pdb = missing_openmm
    snapshot_from_testsystem = missing_openmm
    to_openmm_topology = missing_openmm
    trajectory_from_mdtraj = missing_openmm
    trajectory_to_mdtraj = missing_openmm
    Snapshot = missing_openmm
    MDSnapshot = missing_openmm
else:
    from .engine import OpenMMEngine as Engine
    from .tools import (
        empty_snapshot_from_openmm_topology,
        snapshot_from_pdb,
        snapshot_from_testsystem,
        to_openmm_topology,
        trajectory_from_mdtraj,
        trajectory_to_mdtraj
    )

    from . import features

    from .snapshot import Snapshot, MDSnapshot
    from . import topology

from openpathsampling.engines import NoEngine, SnapshotDescriptor
