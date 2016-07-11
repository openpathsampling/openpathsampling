"""
Attributes
----------
state : :class:`openpathsampling.DynamicsEngine`
    referenec to the engine used to generate the snapshot
"""

variables = ['state']


def netcdfplus_init(store):
    store.create_variable(
        'state',
        'int',
        description="the current discrete state index"
    )
