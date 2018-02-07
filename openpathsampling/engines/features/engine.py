"""
Attributes
----------
engine : :class:`openpathsampling.DynamicsEngine`
    reference to the engine used to generate the snapshot
"""

variables = ['engine']

schema_entries = [('engine', 'uuid')]


def netcdfplus_init(store):
    store.create_variable(
        'engine',
        'obj.engines',
        description="reference to the engine for the current simulation box"
    )


@property
def descriptor(snapshot):
    return snapshot.engine.descriptor


@property
def topology(snapshot):
    return snapshot.engine.topology
