"""
Attributes
----------
engine : :class:`openpathsampling.DynamicsEngine`
    referenec to the engine used to generate the snapshot
"""

attributes = ['engine', 'topology']


def netcdfplus_init(store):
    store.create_variable(
        'engine',
        'obj.engines',
        description="reference to the engine for the current simulation box"
    )


def topology(snapshot):
    return snapshot.engine.topology