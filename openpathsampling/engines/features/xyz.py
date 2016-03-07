attributes = ['xyz']
required = ['coordinates']


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
