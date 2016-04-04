variables = ['box_vectors']
numpy = ['box_vectors']


def netcdfplus_init(store):
    store.create_variable('box_vectors', 'numpy.float32',
                        dimensions=('spatial', 'spatial'),
                        description="box_vectors for the current simulation box",
                        chunksizes=(1, 'atom', 'spatial')
                        )
