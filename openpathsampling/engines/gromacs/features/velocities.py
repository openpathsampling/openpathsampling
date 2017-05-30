@property
def velocities(snapshot):
    """
    """
    if snapshot._velocities is None:
        snapshot.load_details()

    return snapshot._velocities
