attributes = ['xyz']


def xyz(snapshot):
    """
    Coordinates without dimensions.

    SERIOUS PROBLEM: whether .coordinates returns a u.Quantity or jut
    numpy array depending on situation (at least for ToyDynamics). This
    is bad.
    """
    import simtk.unit as u

    coord = snapshot.coordinates
    if type(coord) is u.Quantity:
        return coord._value
    else:
        return coord
