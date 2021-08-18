try:
    import hoomd
except ImportError:
    pass
else:
    from openpathsampling.engines.features import (
        attach_features,
        box_vectors,
        coordinates,
        engine,
        kinetics,
        statics,
        topology,
        velocities,
    )
    from openpathsampling.engines.features.shared import (
        KineticContainer,
        StaticContainer,
    )
