variables = ['box_vectors']
numpy = ['box_vectors']

dimensions = ['n_spatial']


def netcdfplus_init(store):
    store.create_variable(
        'box_vectors', 'numpy.float32',
        dimensions=('n_spatial', 'n_spatial'),
        description="box_vectors for the current simulation box",
        chunksizes=('n_spatial', 'n_spatial'))
