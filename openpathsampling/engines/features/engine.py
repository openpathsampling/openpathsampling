"""
Attributes
----------
engine : :class:`openpathsampling.DynamicsEngine`
    referenec to the engine used to generate the snapshot
"""

variables = ['engine']

def netcdfplus_init(store):
    store.create_variable(
        'engine',
        'obj.engines',
        description="reference to the engine for the current simulation box"
    )


@property
def dim(snapshot):
    return snapshot.engine.dimensions
