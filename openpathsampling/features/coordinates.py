attributes = ['coordinates']


def netcdfplus_init(store):
    store.create_variable('coordinates', 'numpy.float32',
                        dimensions=('atom', 'spatial'),
                        description="coordinate of atom '{ix[1]}' in dimension " +
                                   "'{ix[2]}' of configuration '{ix[0]}'.",
                        chunksizes=(1, 'atom', 'spatial')
                        )
