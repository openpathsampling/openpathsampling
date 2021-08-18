def missing_hoomd(*args, **kwargs):
    raise RuntimeError("Install HOOMD-blue >= 3.0 to use this feature")


try:
    import hoomd
except ImportError:
    HAS_HOOMD = False
    Engine = missing_hoomd
    empty_snapshot_from_openmm_topology = missing_hoomd
    snapshot_from_pdb = missing_hoomd
    snapshot_from_testsystem = missing_hoomd
    to_openmm_topology = missing_hoomd
    trajectory_from_mdtraj = missing_hoomd
    trajectory_to_mdtraj = missing_hoomd
    Snapshot = missing_hoomd
    MDSnapshot = missing_hoomd
else:
    from .engine import HOOMDEngine as Engine

    from . import features

    from .snapshot import Snapshot, MDSnapshot
    from . import topology

from openpathsampling.engines import NoEngine, SnapshotDescriptor
