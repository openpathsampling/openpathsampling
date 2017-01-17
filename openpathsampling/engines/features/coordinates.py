"""
Attributes
----------
coordinates : numpy.ndarray, shape=(atoms, 3), dtype=numpy.float32
    atomic coordinates
"""

variables = ['coordinates']
numpy = ['coordinates']

dimensions = ['n_atoms', 'n_spatial']


def netcdfplus_init(store):
    store.create_variable(
        'coordinates', 'numpy.float32',
        dimensions=('n_atoms', 'n_spatial'),
        description="coordinate of atom '{ix[1]}' in dimension " +
                   "'{ix[2]}' of configuration '{ix[0]}'.",
        chunksizes=('n_atoms', 'n_spatial'))


@property
def xyz(snapshot):
    """
    Returns
    -------
    xyz : numpy.ndarray, shape=(atoms, 3), dtype=numpy.float32
        atomic coordinates without dimensions. Be careful.

    """
    import simtk.unit as u

    coord = snapshot.coordinates
    if type(coord) is u.Quantity:
        return coord._value
    else:
        return coord
