default_none = ['_xyz']

@property
def xyz(snapshot):
    """
    Returns
    -------
    xyz : numpy.ndarray, shape=(atoms, 3), dtype=numpy.float32
        atomic coordinates without dimensions.
    """
    if snapshot._xyz is None:
        snapshot.load_details()

    return snapshot._xyz

@property
def coordinates(snapshot):
    """
    """
    return snapshot.xyz
