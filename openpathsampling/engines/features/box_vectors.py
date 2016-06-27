variables = ['box_vectors']
numpy = ['box_vectors']

dimensions = ['spatial']


def netcdfplus_init(store):
    store.create_variable('box_vectors', 'numpy.float32',
                        dimensions=('spatial', 'spatial'),
                        description="box_vectors for the current simulation box",
                        chunksizes=(1, 'spatial', 'spatial')
                        )
