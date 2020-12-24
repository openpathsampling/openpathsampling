try:
    import simtk.openmm
    import simtk.openmm.app
except ImportError:
    pass
else:
    from openpathsampling.engines.features import *
    from openpathsampling.engines.features.shared import StaticContainer, KineticContainer
    from . import masses
    from . import instantaneous_temperature
    from . import traj_quantities
