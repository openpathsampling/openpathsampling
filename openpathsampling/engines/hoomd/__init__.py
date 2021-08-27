def missing_hoomd(*args, **kwargs):
    raise RuntimeError("Install HOOMD-blue >= 3.0.0 to use this feature")


try:
    import hoomd
    if hoomd.version.version < "3.0.0":
        raise ValueError("HOOMD-blue >= 3.0.0 is required.")
except (ImportError, AttributeError, ValueError):
    HAS_HOOMD = False
    Engine = missing_hoomd
    Snapshot = missing_hoomd
else:
    from .engine import HOOMDEngine as Engine
    from . import features
    from .snapshot import Snapshot
