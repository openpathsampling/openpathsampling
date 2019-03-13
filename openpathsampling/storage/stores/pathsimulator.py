from openpathsampling.pathsimulators import PathSimulator
from openpathsampling.netcdfplus import NamedObjectStore


class PathSimulatorStore(NamedObjectStore):
    def __init__(self):
        super(PathSimulatorStore, self).__init__(
            PathSimulator,
        )

    def to_dict(self):
        return {}

    def _load(self, idx):
        obj = super(PathSimulatorStore, self)._load(idx)
        if hasattr(obj, 'storage'):
            obj.storage = self.storage

        return obj
