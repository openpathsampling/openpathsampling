from openpathsampling.netcdfplus import AttributeStore


class CVStore(AttributeStore):
    """
    ObjectStore to store a dict with StorableObject : value
    """
    def __init__(self):
        super(CVStore, self).__init__()
