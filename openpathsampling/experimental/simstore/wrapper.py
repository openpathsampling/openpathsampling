from openpathsampling.netcdfplus import StorableNamedObject

class SimStoreWrapper(StorableNamedObject):
    """Wrapper to add UUID to and JSON-serializable content.

    Parameters
    ----------
    content : Any
        the JSON-serializable thing to wrap
    """
    def __init__(self, content):
        super().__init__()
        self.content = content
