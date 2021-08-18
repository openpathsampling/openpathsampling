import openpathsampling.engines.topology as topology
from openpathsampling.deprecations import OPENMM_MDTRAJTOPOLOGY

# Mimicked from collections to warn on import
def __getattr__(name):
    if name == "MDTrajTopology":
        # Needs to be stacklevel=3 to surface
        OPENMM_MDTRAJTOPOLOGY.warn(stacklevel=3)
        obj = getattr(topology, name)
        globals()[name] = obj
        return obj
    raise AttributeError(
        "module '{}' has no attribute '{}'".format(__name__, name)
    )
