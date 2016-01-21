_variables = ['box_vectors']

attributes = ['box_vectors']
attributes_minus = []
attributes_not = []

properties = []

def _init(store):
    store.create_variable('box_vectors', 'numpy.float32',
                        dimensions=('spatial', 'spatial'),
                        description="box_vectors for the current simulation box",
                        chunksizes=(1, 'atom', 'spatial')
                        )
