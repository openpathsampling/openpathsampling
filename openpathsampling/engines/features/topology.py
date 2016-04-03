"""
Attributes
----------
topology : :class:`openpathsampling.Topolgy`
    object containing a description of the snapshots topology
"""


variables = ['topology']


def netcdfplus_init(store):
    store.create_variable('topology', 'obj.topologies',
                        description="reference to the topology for the current simulation box"
                        )
