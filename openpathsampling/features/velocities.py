"""
velocities : simtk.unit.Quantity wrapping Nx3 np array of dimension length
    atomic velocities (default: None)
"""

attributes = ['velocities']
minus = ['velocities']


def netcdfplus_init(store):

    store.create_variable('velocities', 'numpy.float32',
                        dimensions=('atom', 'spatial'),
                        description="the velocity of atom 'atom' in dimension " +
                                   "'coordinate' of momentum 'momentum'.",
                        chunksizes=(1, 'atom', 'spatial')
                        )