"""
Attributes
----------
topology : :class:`openpathsampling.Topology`
    object containing a description of the snapshot's topology
"""


variables = ['topology']


def netcdfplus_init(store):
    store.create_variable('topology', 'obj.topologies',
                        description="reference to the topology for the current simulation box"
                        )
