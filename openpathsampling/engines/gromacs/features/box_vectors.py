@property
def box_vectors(snapshot):
    """
    """
    if snapshot._box_vectors is None:
        snapshot.load_details()

    return snapshot._box_vectors
