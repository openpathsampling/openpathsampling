"""
velocities : simtk.unit.Quantity wrapping Nx3 np array of dimension length
    atomic velocities (default: None)
"""

_variables = ['velocities']

attributes = ['velocities']
attributes_minus = ['velocities']
attributes_not = []

properties = []

def _init(store):

    store.create_variable('velocities', 'numpy.float32',
                        dimensions=('atom', 'spatial'),
                        description="the velocity of atom 'atom' in dimension " +
                                   "'coordinate' of momentum 'momentum'.",
                        chunksizes=(1, 'atom', 'spatial')
                        )