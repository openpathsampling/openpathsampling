def requires_mdtraj(*args, **kwargs):  # pragma: no cover
    raise RuntimeError("This requires MDTraj, which is not installed")

from openpathsampling.integration_tools import HAS_MDTRAJ

if not HAS_MDTRAJ:
    Engine = requires_mdtraj
    ExternalMDSnapshot = requires_mdtraj
    snapshot_from_gro = requires_mdtraj
else:
    from .engine import GromacsEngine as Engine
    from .engine import ExternalMDSnapshot
    from .engine import snapshot_from_gro

