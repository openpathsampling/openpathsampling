try:
    import hoomd
except ImportError:
    pass
else:
    from openpathsampling.engines.features import statics, coordinates, velocities, kinetics, box_vectors, topology, engine
    from openpathsampling.engines.features.shared import StaticContainer, KineticContainer
    from . import masses
    from . import instantaneous_temperature
    from . import traj_quantities
