_variables = ['velocities']


def _init(store):

    store.create_variable('velocities', 'numpy.float32',
                        dimensions=('atom', 'spatial'),
                        description="the velocity of atom 'atom' in dimension " +
                                   "'coordinate' of momentum 'momentum'.",
                        chunksizes=(1, 'atom', 'spatial')
                        )