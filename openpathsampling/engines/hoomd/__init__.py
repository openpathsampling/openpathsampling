def missing_hoomd(*args, **kwargs):
    raise RuntimeError("Install HOOMD-blue >= 3.0 to use this feature")


try:
    import hoomd
except ImportError:
    HAS_HOOMD = False
    Engine = missing_hoomd
    Snapshot = missing_hoomd
else:
    from .engine import HOOMDEngine as Engine
    from . import features
    from .snapshot import Snapshot
